"""
plot_radial_profiles.py - Inciso 1.4
Generates: Radial profiles for fresh inward-moving particles.

Per the enunciado:
- S = distance from center of fixed particle (= origin)
- dS = 0.2 m concentric shell width
- Select only fresh particles Pf^in with Rj · vj < 0 (velocity toward center)
- Compute:
  - <rho_f^in>(S): density = count / shell_area, averaged over saved snapshots
  - <vf_in>(S): radial velocity vf_in_j = (Rj · vj) / |Rj|, averaged
    by pooling all saved particles in each shell before averaging
  - J_in(S) = <rho_f^in>(S) * |<vf_in>(S)|
- Plot 3 curves: <rho_f^in>(S), |<vf_in>(S)|, J_in(S)
- Near obstacle S in [2,3]: plot J_in, <rho_f^in>, <vf_in> as functions of N

Usage:
    python graphics/plot_radial_profiles.py data/sim_300N_file.txt
    python graphics/plot_radial_profiles.py data/
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analysis_cache import group_entries_by_N, load_analysis_entries, load_analysis_file
from plot_fraction_used import MANUAL_T_EST_BY_N
from plot_style import apply_plot_style, get_distinct_series_styles, format_y_axis


DISPLAY_PROFILE_N_VALUES = {100, 300, 500, 800}
NEAR_OBSTACLE_S_RANGE = (2.0, 3.0)
ZOOM_S_RANGE = (2.0, 5.0)
STATIONARY_ACCUMULATORS_CACHE = {}


apply_plot_style()


def std_to_sem(std_values, sample_counts):
    """Convert standard deviation arrays to standard error arrays (std/sqrt(n))."""
    std_values = np.asarray(std_values, dtype=float)
    sample_counts = np.asarray(sample_counts, dtype=float)
    sem = np.zeros_like(std_values, dtype=float)
    valid = sample_counts > 0.0
    sem[valid] = std_values[valid] / np.sqrt(sample_counts[valid])
    return sem


def propagate_flux_error(rho_mean, v_mean, rho_err, v_err):
    """First-order uncertainty propagation for J = rho * |v|."""
    rho_mean = np.asarray(rho_mean, dtype=float)
    v_mean = np.asarray(v_mean, dtype=float)
    rho_err = np.asarray(rho_err, dtype=float)
    v_err = np.asarray(v_err, dtype=float)
    return np.sqrt((np.abs(v_mean) * rho_err) ** 2 + (rho_mean * v_err) ** 2)


def reduce_shell_band(values, errors, shell_mask):
    """Average one observable over a shell band and combine independent errors."""
    band_values = np.asarray(values, dtype=float)[shell_mask]
    band_errors = np.asarray(errors, dtype=float)[shell_mask]
    if band_values.size == 0:
        return np.nan, np.nan
    mean_value = float(np.mean(band_values))
    combined_error = float(np.sqrt(np.sum(band_errors ** 2)) / band_values.size)
    return mean_value, combined_error


def pick_shell_band_mask(S_centers, s_min, s_max):
    """Pick shells inside [s_min, s_max], fallback to the first shell with center >= s_min."""
    S_centers = np.asarray(S_centers, dtype=float)
    mask = (S_centers >= s_min) & (S_centers <= s_max)
    if np.any(mask):
        return mask
    fallback_idx = np.where(S_centers >= s_min)[0]
    if fallback_idx.size == 0:
        fallback_idx = np.array([len(S_centers) - 1])
    mask = np.zeros_like(S_centers, dtype=bool)
    mask[fallback_idx[0]] = True
    return mask


def parse_metadata_tokens(line, metadata):
    """Parse numeric key=value tokens from one header line."""
    stripped = line.strip().lstrip('# ')
    for part in stripped.split():
        if '=' not in part:
            continue
        key, val = part.split('=', 1)
        try:
            metadata[key] = float(val)
        except ValueError:
            metadata[key] = val


def compute_stationary_accumulators_from_file(filepath, t_est, dS=0.2):
    """
    Stream one sim file and accumulate stationary radial statistics.

    For each shell we keep:
    - one density sample per snapshot: the inward fresh-particle count in that shell
    - one velocity sample per particle: each individual inward radial velocity
    """
    metadata = {}
    shell_edges = None
    S_centers = None
    shell_areas = None
    count_sum = None
    count_sq_sum = None
    velocity_sample_count = None
    velocity_sum = None
    velocity_sq_sum = None
    snapshot_count = 0.0
    N = None

    with open(filepath, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            trimmed = line.strip()
            if not trimmed:
                continue

            if trimmed.startswith('#'):
                parse_metadata_tokens(trimmed, metadata)
                continue

            if not trimmed.startswith('S '):
                continue

            if N is None:
                N = int(metadata.get('N', 0))
                if N <= 0:
                    raise ValueError(f"Invalid or missing N in metadata for {filepath}")
                shell_edges, S_centers, shell_areas = build_shell_geometry(metadata, dS)
                count_sum = np.zeros(len(S_centers), dtype=float)
                count_sq_sum = np.zeros(len(S_centers), dtype=float)
                velocity_sample_count = np.zeros(len(S_centers), dtype=float)
                velocity_sum = np.zeros(len(S_centers), dtype=float)
                velocity_sq_sum = np.zeros(len(S_centers), dtype=float)

            t_snapshot = float(trimmed.split()[1])
            include_snapshot = (t_est is None) or (t_snapshot >= t_est)
            snapshot_counts = np.zeros(len(S_centers), dtype=float) if include_snapshot else None

            for _ in range(N):
                particle_line = f.readline()
                if not particle_line:
                    raise EOFError(f"Unexpected EOF while reading snapshot at t={t_snapshot} in {filepath}")

                if not include_snapshot:
                    continue

                parts = particle_line.strip().split()
                if len(parts) < 6 or parts[5] != 'F':
                    continue

                x = float(parts[1])
                y = float(parts[2])
                vx = float(parts[3])
                vy = float(parts[4])
                S = np.sqrt(x * x + y * y)
                rdotv = x * vx + y * vy
                if rdotv >= 0:
                    continue
                if S < shell_edges[0] - 1e-12 or S > shell_edges[-1] + 1e-12:
                    continue

                shell_idx = np.searchsorted(shell_edges, S, side='right') - 1
                shell_idx = min(max(shell_idx, 0), len(S_centers) - 1)
                vf_in = rdotv / S if S > 1e-10 else 0.0
                snapshot_counts[shell_idx] += 1.0
                count_sum[shell_idx] += 1.0
                velocity_sample_count[shell_idx] += 1.0
                velocity_sum[shell_idx] += vf_in
                velocity_sq_sum[shell_idx] += vf_in * vf_in

            if include_snapshot:
                count_sq_sum += snapshot_counts * snapshot_counts
                snapshot_count += 1.0

    if S_centers is None:
        raise ValueError(f"No snapshots found in {filepath}")

    return (
        S_centers,
        shell_areas,
        snapshot_count,
        count_sum,
        count_sq_sum,
        velocity_sample_count,
        velocity_sum,
        velocity_sq_sum,
    )


def parse_simulation_file(filepath):
    """Parse simulation file. Returns metadata, snapshots, events."""
    metadata = {}
    snapshots = []
    events = []
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    i = 0
    total = len(lines)
    
    while i < total and lines[i].startswith('#'):
        line = lines[i].strip().lstrip('# ')
        if 'N=' in line:
            parts = line.split()
            for part in parts:
                if '=' in part:
                    key, val = part.split('=')
                    try:
                        metadata[key] = float(val)
                    except ValueError:
                        metadata[key] = val
        i += 1
    
    N = int(metadata.get('N', 0))
    
    while i < total:
        line = lines[i].strip()
        if line.startswith('S '):
            time = float(line.split()[1])
            particles = []
            for j in range(N):
                i += 1
                if i >= total:
                    break
                parts = lines[i].strip().split()
                if len(parts) >= 6:
                    particles.append({
                        'id': int(parts[0]),
                        'x': float(parts[1]),
                        'y': float(parts[2]),
                        'vx': float(parts[3]),
                        'vy': float(parts[4]),
                        'state': parts[5]
                    })
            snapshots.append((time, particles))
        elif line.startswith('E '):
            parts = line.split()
            events.append((float(parts[1]), int(parts[2])))
        i += 1
    
    return metadata, snapshots, events


def build_shell_geometry(metadata, dS=0.2):
    """Return physically valid shell edges, centers, and areas."""
    R_enclosure = float(metadata['R_enclosure'])
    r0 = float(metadata['r0'])
    r_particle = float(metadata['r'])

    # Particle centers can only move between these radii.
    S_min = r0 + r_particle
    S_max = R_enclosure - r_particle

    shell_edges = np.arange(S_min, S_max + 1e-9, dS)
    if shell_edges[-1] < S_max:
        shell_edges = np.append(shell_edges, S_max)

    S_centers = (shell_edges[:-1] + shell_edges[1:]) / 2.0
    shell_areas = np.pi * (shell_edges[1:] ** 2 - shell_edges[:-1] ** 2)
    return shell_edges, S_centers, shell_areas


def compute_radial_profiles(snapshots, metadata, dS=0.2):
    """
    Compute radial profiles for fresh inward-moving particles.

    This helper mirrors the cache-based logic: it uses only the saved
    snapshots and pools all particles observed in each shell before averaging.
    """
    S_centers, shell_areas, snapshot_count, count_sum, velocity_sum = \
        compute_radial_accumulators(snapshots, metadata, dS)
    rho_mean, v_mean, J_in = derive_profiles_from_accumulators(
        snapshot_count,
        shell_areas,
        count_sum,
        velocity_sum,
    )
    return S_centers, rho_mean, v_mean, J_in


def compute_radial_accumulators(snapshots, metadata, dS=0.2):
    """Accumulate shell counts/velocity sums from a snapshot list."""
    shell_edges, S_centers, shell_areas = build_shell_geometry(metadata, dS)
    n_shells = len(S_centers)
    count_sum = np.zeros(n_shells)
    velocity_sum = np.zeros(n_shells)

    for _, particles in snapshots:
        for p in particles:
            if p['state'] != 'F':
                continue

            x, y = p['x'], p['y']
            vx, vy = p['vx'], p['vy']
            S = np.sqrt(x * x + y * y)

            rdotv = x * vx + y * vy
            if rdotv >= 0:
                continue

            if S < shell_edges[0] - 1e-12 or S > shell_edges[-1] + 1e-12:
                continue

            shell_idx = np.searchsorted(shell_edges, S, side='right') - 1
            shell_idx = min(max(shell_idx, 0), n_shells - 1)

            vf_in = rdotv / S if S > 1e-10 else 0.0
            count_sum[shell_idx] += 1
            velocity_sum[shell_idx] += vf_in

    return S_centers, shell_areas, float(len(snapshots)), count_sum, velocity_sum


def filter_stationary_snapshots(snapshots, t_est):
    """Keep only snapshots in the stationary regime t >= t_est."""
    if t_est is None:
        return snapshots
    return [(t, particles) for t, particles in snapshots if t >= t_est]


def extract_stationary_accumulators(entry, dS=0.2):
    """Build radial accumulators for one realization using only t >= t_est snapshots."""
    N = int(entry["metadata"]["N"])
    t_est = MANUAL_T_EST_BY_N.get(N)
    key = (
        entry["source_path"],
        float(dS),
        t_est,
        entry.get("source_mtime_ms"),
        entry.get("source_size"),
    )
    if key not in STATIONARY_ACCUMULATORS_CACHE:
        STATIONARY_ACCUMULATORS_CACHE[key] = compute_stationary_accumulators_from_file(
            entry["source_path"],
            t_est,
            dS,
        )
    return STATIONARY_ACCUMULATORS_CACHE[key], t_est


def aggregate_radial_profiles(entries, dS=0.2):
    """
    Aggregate radial profiles over multiple realizations.

    Density stats use one sample per snapshot and velocity stats use one sample
    per particle, both pooled across every stationary snapshot in every
    realization.
    """
    stationary_data = [extract_stationary_accumulators(entry, dS) for entry in entries]
    accumulators = [item[0] for item in stationary_data]
    S_centers = accumulators[0][0]
    shell_areas = accumulators[0][1]

    pooled_snapshot_count = sum(acc[2] for acc in accumulators)
    pooled_count_sum = np.sum([acc[3] for acc in accumulators], axis=0)
    pooled_count_sq_sum = np.sum([acc[4] for acc in accumulators], axis=0)
    pooled_velocity_sample_count = np.sum([acc[5] for acc in accumulators], axis=0)
    pooled_vf_in_sum = np.sum([acc[6] for acc in accumulators], axis=0)
    pooled_vf_in_sq_sum = np.sum([acc[7] for acc in accumulators], axis=0)
    rho_mean, rho_std, v_mean, v_std, J_mean = derive_profiles_from_stationary_samples(
        pooled_snapshot_count,
        shell_areas,
        pooled_count_sum,
        pooled_count_sq_sum,
        pooled_velocity_sample_count,
        pooled_vf_in_sum,
        pooled_vf_in_sq_sum,
    )

    rho_err = std_to_sem(rho_std, np.full_like(rho_std, pooled_snapshot_count, dtype=float))
    v_err = std_to_sem(v_std, pooled_velocity_sample_count)
    J_err = propagate_flux_error(rho_mean, v_mean, rho_err, v_err)

    return (
        S_centers,
        rho_mean,
        rho_err,
        v_mean,
        v_err,
        J_mean,
        J_err,
    )


def extract_profile_statistics_for_entries(entries, dS=0.2):
    """Convenience wrapper: pooled profile means + SEM errors for one N."""
    S_centers, rho_mean, rho_err, v_mean, v_err, J_mean, J_err = aggregate_radial_profiles(entries, dS)
    return {
        "S_centers": S_centers,
        "rho_mean": rho_mean,
        "rho_err": rho_err,
        "v_mean": v_mean,
        "v_err": v_err,
        "J_mean": J_mean,
        "J_err": J_err,
    }


def collect_profiles_by_N(files_by_N, dS=0.2):
    """Compute pooled radial statistics for each N once."""
    profile_by_N = {}
    for N in sorted(files_by_N.keys()):
        profile_by_N[N] = extract_profile_statistics_for_entries(files_by_N[N], dS)
    return profile_by_N


def derive_profiles_from_accumulators(snapshot_count, shell_areas, count_sum, vf_in_sum):
    """Build rho, v and J from pooled shell accumulators."""
    count_sum = np.array(count_sum, dtype=float)
    vf_in_sum = np.array(vf_in_sum, dtype=float)
    shell_areas = np.array(shell_areas, dtype=float)

    rho = np.zeros_like(count_sum, dtype=float)
    v = np.full_like(count_sum, np.nan, dtype=float)
    J = np.zeros_like(count_sum, dtype=float)

    if snapshot_count > 0:
        rho = count_sum / (float(snapshot_count) * shell_areas)

    nonzero = count_sum > 1e-12
    v[nonzero] = vf_in_sum[nonzero] / count_sum[nonzero]
    J[nonzero] = rho[nonzero] * np.abs(v[nonzero])
    return rho, v, J


def derive_profiles_from_stationary_samples(
    snapshot_count,
    shell_areas,
    count_sum,
    count_sq_sum,
    velocity_sample_count,
    velocity_sum,
    velocity_sq_sum,
):
    """
    Build pooled profile means/stds from stationary samples.

    Density uses one sample per snapshot: the shell count divided by shell area.
    Velocity uses one sample per particle: each inward radial velocity value.
    """
    count_sum = np.array(count_sum, dtype=float)
    count_sq_sum = np.array(count_sq_sum, dtype=float)
    velocity_sample_count = np.array(velocity_sample_count, dtype=float)
    velocity_sum = np.array(velocity_sum, dtype=float)
    velocity_sq_sum = np.array(velocity_sq_sum, dtype=float)
    shell_areas = np.array(shell_areas, dtype=float)

    rho_mean = np.zeros_like(count_sum, dtype=float)
    rho_std = np.zeros_like(count_sum, dtype=float)
    if snapshot_count > 0:
        mean_count = count_sum / float(snapshot_count)
        mean_count_sq = count_sq_sum / float(snapshot_count)
        var_count = np.maximum(mean_count_sq - mean_count * mean_count, 0.0)
        rho_mean = mean_count / shell_areas
        rho_std = np.sqrt(var_count) / shell_areas

    v_mean = np.full_like(count_sum, np.nan, dtype=float)
    v_std = np.full_like(count_sum, np.nan, dtype=float)
    nonzero_velocity = velocity_sample_count > 1e-12
    if np.any(nonzero_velocity):
        mean_velocity = velocity_sum[nonzero_velocity] / velocity_sample_count[nonzero_velocity]
        mean_velocity_sq = velocity_sq_sum[nonzero_velocity] / velocity_sample_count[nonzero_velocity]
        var_velocity = np.maximum(mean_velocity_sq - mean_velocity * mean_velocity, 0.0)
        v_mean[nonzero_velocity] = mean_velocity
        v_std[nonzero_velocity] = np.sqrt(var_velocity)

    J_mean = np.zeros_like(count_sum, dtype=float)
    J_mean[nonzero_velocity] = rho_mean[nonzero_velocity] * np.abs(v_mean[nonzero_velocity])
    return rho_mean, rho_std, v_mean, v_std, J_mean


def extract_radial_profile_accumulators(entry):
    """Return cached shell accumulators for one realization."""
    radial = entry["radial_profiles"]
    return (
        np.array(radial["S_centers"], dtype=float),
        np.array(radial["shell_areas"], dtype=float),
        float(radial["snapshot_count"]),
        np.array(radial["count_sum"], dtype=float),
        np.array(radial["vf_in_sum"], dtype=float),
    )


def extract_radial_profiles(entry):
    """Return cached radial profile arrays for one realization."""
    S_centers, shell_areas, snapshot_count, count_sum, vf_in_sum = \
        extract_radial_profile_accumulators(entry)
    rho, v, J = derive_profiles_from_accumulators(
        snapshot_count,
        shell_areas,
        count_sum,
        vf_in_sum,
    )
    return S_centers, rho, v, J


def save_profile_plot(
    x_values,
    y_values,
    outpath,
    xlabel,
    ylabel,
    *,
    fmt,
    color,
    markerfacecolor,
    markeredgecolor,
    markersize,
    linewidth,
    scientific_y=False,
    yerr=None,
    capsize=0,
    xticks=None,
):
    """Save one radial-profile-related plot, optionally with vertical error bars."""
    fig, ax = plt.subplots(figsize=(10, 7))

    if yerr is None:
        ax.plot(
            x_values,
            y_values,
            fmt,
            color=color,
            markerfacecolor=markerfacecolor,
            markeredgecolor=markeredgecolor,
            markersize=markersize,
            linewidth=linewidth,
        )
    else:
        ax.errorbar(
            x_values,
            y_values,
            yerr=yerr,
            fmt=fmt,
            capsize=capsize,
            color=color,
            markerfacecolor=markerfacecolor,
            markeredgecolor=markeredgecolor,
            markersize=markersize,
            linewidth=linewidth,
        )

    ax.set_xlabel(xlabel, fontsize=17)
    ax.set_ylabel(ylabel, fontsize=17)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', labelsize=14)

    if xticks is not None:
        ax.set_xticks(xticks)

    format_y_axis(ax, yerr=yerr)

    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()


def save_radial_profile_figures(
    N,
    S_centers,
    rho_values,
    v_values,
    J_values,
    output_dir,
    *,
    rho_err=None,
    v_err=None,
    J_err=None,
):
    """Save one file per radial observable for a given N."""
    outputs = []

    rho_outpath = os.path.join(output_dir, f"inciso_1_4_rho_profile_N{N}.png")
    save_profile_plot(
        S_centers,
        rho_values,
        rho_outpath,
        'Distancia radial $S$ [m]',
        r'$\langle \rho_f^{in} \rangle$ [partículas/m²]',
        fmt='o-',
        color='#1B5E20',
        markerfacecolor='#4CAF50',
        markeredgecolor='#1B5E20',
        markersize=4,
        linewidth=1.4,
        scientific_y=True,
        yerr=rho_err,
        capsize=3,
    )
    outputs.append(rho_outpath)

    v_outpath = os.path.join(output_dir, f"inciso_1_4_velocity_profile_N{N}.png")
    save_profile_plot(
        S_centers,
        np.abs(v_values),
        v_outpath,
        'Distancia radial $S$ [m]',
        r'$|\langle v_f^{in} \rangle|$ [m/s]',
        fmt='s-',
        color='#0D47A1',
        markerfacecolor='#2196F3',
        markeredgecolor='#0D47A1',
        markersize=4,
        linewidth=1.4,
        yerr=v_err,
        capsize=3,
    )
    outputs.append(v_outpath)

    j_outpath = os.path.join(output_dir, f"inciso_1_4_flux_profile_N{N}.png")
    save_profile_plot(
        S_centers,
        J_values,
        j_outpath,
        'Distancia radial $S$ [m]',
        r'$J_{in}$ [particles/(m·s)]',
        fmt='D-',
        color='#BF360C',
        markerfacecolor='#FF5722',
        markeredgecolor='#BF360C',
        markersize=4,
        linewidth=1.4,
        scientific_y=True,
        yerr=J_err,
        capsize=3,
    )
    outputs.append(j_outpath)

    return outputs


def save_s2_vs_n_figures(
    N_values,
    rho_values,
    v_values,
    J_values,
    output_dir,
    *,
    rho_err,
    v_err,
    J_err,
):
    """Save one file per observable averaged in the near-obstacle band S in [2,3]."""
    outputs = []

    rho_outpath = os.path.join(output_dir, "inciso_1_4_S2_S3_rho_vs_N.png")
    save_profile_plot(
        N_values,
        rho_values,
        rho_outpath,
        'Número de partículas $N$',
        r'$\langle \rho_f^{in} \rangle$ [part/m²] en $S \in [2,3]$ m',
        fmt='o-',
        color='#1B5E20',
        markerfacecolor='#4CAF50',
        markeredgecolor='#1B5E20',
        markersize=8,
        linewidth=2.0,
        scientific_y=True,
        yerr=rho_err,
        capsize=5,
        xticks=N_values,
    )
    outputs.append(rho_outpath)

    v_outpath = os.path.join(output_dir, "inciso_1_4_S2_S3_velocity_vs_N.png")
    save_profile_plot(
        N_values,
        v_values,
        v_outpath,
        'Número de partículas $N$',
        r'$|\langle v_f^{in} \rangle|$ [m/s] en $S \in [2,3]$ m',
        fmt='s-',
        color='#0D47A1',
        markerfacecolor='#2196F3',
        markeredgecolor='#0D47A1',
        markersize=8,
        linewidth=2.0,
        yerr=v_err,
        capsize=5,
        xticks=N_values,
    )
    outputs.append(v_outpath)

    j_outpath = os.path.join(output_dir, "inciso_1_4_S2_S3_flux_vs_N.png")
    save_profile_plot(
        N_values,
        J_values,
        j_outpath,
        'Número de partículas $N$',
        r'$J_{in}$ [particles/(m·s)] en $S \in [2,3]$ m',
        fmt='D-',
        color='#BF360C',
        markerfacecolor='#FF5722',
        markeredgecolor='#BF360C',
        markersize=8,
        linewidth=2.0,
        scientific_y=True,
        yerr=J_err,
        capsize=5,
        xticks=N_values,
    )
    outputs.append(j_outpath)

    return outputs


def compute_zoom_y_limits(profile_by_N, key_mean, key_err, x_min, x_max):
    """Compute y-limits from samples inside the zoomed x-window."""
    y_min = np.inf
    y_max = -np.inf

    for stats in profile_by_N.values():
        S_centers = np.asarray(stats["S_centers"], dtype=float)
        mask = (S_centers >= x_min) & (S_centers <= x_max)
        if not np.any(mask):
            continue
        y = np.asarray(stats[key_mean], dtype=float)
        if key_mean == "v_mean":
            y = np.abs(y)
        y_err = np.asarray(stats[key_err], dtype=float)
        low = np.nanmin(y[mask] - y_err[mask])
        high = np.nanmax(y[mask] + y_err[mask])
        y_min = min(y_min, float(low))
        y_max = max(y_max, float(high))

    if not np.isfinite(y_min) or not np.isfinite(y_max):
        return None

    if y_max - y_min < 1e-12:
        margin = 0.1 * max(abs(y_max), 1.0)
    else:
        margin = 0.08 * (y_max - y_min)

    lower = max(0.0, y_min - margin)
    upper = y_max + margin
    if upper <= lower:
        upper = lower + 1.0
    return lower, upper


def save_profiles_by_n_figures(profile_by_N, output_dir, zoom_range=None):
    """Save rho, |v| and J profiles in the same figure for multiple N."""
    os.makedirs(output_dir, exist_ok=True)
    N_values = sorted(profile_by_N.keys())
    zoom_suffix = "" if zoom_range is None else "_zoom_S2_S5"
    outpath = os.path.join(output_dir, f"inciso_1_4_profiles_vs_S_by_N{zoom_suffix}.png")

    if zoom_range is not None:
        fig, ax = plt.subplots(figsize=(10, 7))
        cmap = plt.get_cmap('viridis')
        n_min = float(min(N_values))
        n_max = float(max(N_values))
        if np.isclose(n_min, n_max):
            norm = matplotlib.colors.Normalize(vmin=n_min - 0.5, vmax=n_max + 0.5)
        else:
            norm = matplotlib.colors.Normalize(vmin=n_min, vmax=n_max)

        y_limits = compute_zoom_y_limits(
            profile_by_N,
            "J_mean",
            "J_err",
            zoom_range[0],
            zoom_range[1],
        )

        all_j_err = []
        for N in N_values:
            stats = profile_by_N[N]
            S_centers = stats["S_centers"]
            all_j_err.append(np.asarray(stats["J_err"], dtype=float))
            ax.errorbar(
                S_centers,
                stats["J_mean"],
                yerr=stats["J_err"],
                fmt='o-',
                linewidth=1.8,
                markersize=4.2,
                capsize=2,
                alpha=0.95,
                color=cmap(norm(float(N))),
            )

        ax.set_xlabel('Distancia radial $S$ [m]', fontsize=16)
        ax.set_ylabel(r'$J_{in}$ [particles/(m·s)]', color='black', fontsize=16)
        ax.grid(True, alpha=0.25, linestyle='--')
        ax.tick_params(axis='both', labelsize=13)
        ax.set_xlim(zoom_range[0], zoom_range[1])
        if y_limits is not None:
            ax.set_ylim(y_limits[0], y_limits[1])
        format_y_axis(ax, yerr=np.concatenate(all_j_err) if all_j_err else None)

        scalar_mappable = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
        scalar_mappable.set_array([])
        colorbar = fig.colorbar(scalar_mappable, ax=ax, pad=0.02)
        colorbar.set_label('Numero de particulas $N$', fontsize=13)
        colorbar.ax.tick_params(labelsize=12)

        fig.tight_layout()
        fig.savefig(outpath, dpi=150, bbox_inches='tight')
        plt.close(fig)
        return outpath

    styles = get_distinct_series_styles(len(N_values))

    fig, axes = plt.subplots(1, 3, figsize=(20, 7), sharex=True)
    labels = [
        (r'$\langle \rho_f^{in} \rangle$ [part/m$^2$]', "rho_mean", "rho_err", '#1B5E20'),
        (r'$|\langle v_f^{in} \rangle|$ [m/s]', "v_mean", "v_err", '#0D47A1'),
        (r'$J_{in}$ [particles/(m·s)]', "J_mean", "J_err", '#BF360C'),
    ]

    zoom_y_limits = [None, None, None]
    if zoom_range is not None:
        for idx, (_, key_mean, key_err, _) in enumerate(labels):
            zoom_y_limits[idx] = compute_zoom_y_limits(
                profile_by_N,
                key_mean,
                key_err,
                zoom_range[0],
                zoom_range[1],
            )

    for style, N in zip(styles, N_values):
        stats = profile_by_N[N]
        S_centers = stats["S_centers"]
        for axis, (_, key_mean, key_err, _) in zip(axes, labels):
            y = np.abs(stats[key_mean]) if key_mean == "v_mean" else stats[key_mean]
            axis.errorbar(
                S_centers,
                y,
                yerr=stats[key_err],
                fmt=f'{style["marker"]}{style["linestyle"]}',
                linewidth=1.8,
                markersize=4.5,
                capsize=2,
                alpha=0.95,
                color=style["color"],
                label=f'N={N}',
            )

    yerr_by_axis = [[], [], []]
    for N in N_values:
        stats = profile_by_N[N]
        yerr_by_axis[0].append(np.asarray(stats["rho_err"], dtype=float))
        yerr_by_axis[1].append(np.asarray(stats["v_err"], dtype=float))
        yerr_by_axis[2].append(np.asarray(stats["J_err"], dtype=float))

    for idx, (axis, (ylabel, _, _, _axis_color)) in enumerate(zip(axes, labels)):
        axis.set_xlabel('Distancia radial $S$ [m]', fontsize=16)
        axis.set_ylabel(ylabel, color='black', fontsize=16)
        axis.grid(True, alpha=0.25, linestyle='--')
        axis.tick_params(axis='both', labelsize=13)
        if zoom_range is not None:
            axis.set_xlim(zoom_range[0], zoom_range[1])
            y_limits = zoom_y_limits[idx]
            if y_limits is not None:
                axis.set_ylim(y_limits[0], y_limits[1])
        merged_err = np.concatenate(yerr_by_axis[idx]) if yerr_by_axis[idx] else None
        format_y_axis(axis, yerr=merged_err)

    handles, legends = axes[0].get_legend_handles_labels()
    fig.legend(handles, legends, loc='upper center', ncol=min(len(N_values), 6), fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close(fig)
    return outpath


def save_near_obstacle_vs_scanning_figure(N_values, J_fresh, J_fresh_err, J_scan, J_scan_err, output_dir):
    """Compare J_fresh^in near obstacle with scanning-rate J(N)."""
    os.makedirs(output_dir, exist_ok=True)
    outpath = os.path.join(output_dir, 'inciso_1_4_jfresh_vs_scanning_rate.png')

    fig, ax_left = plt.subplots(figsize=(10, 7))
    ax_right = ax_left.twinx()

    left_line = ax_left.errorbar(
        N_values,
        J_fresh,
        yerr=J_fresh_err,
        fmt='D-',
        capsize=5,
        linewidth=2,
        color='#BF360C',
        markerfacecolor='#FF5722',
        markeredgecolor='#BF360C',
        label=r'$J_{fresh}^{in}(S\in[2,3])$',
    )
    right_line = ax_right.errorbar(
        N_values,
        J_scan,
        yerr=J_scan_err,
        fmt='s--',
        capsize=5,
        linewidth=2,
        color='#6A1B9A',
        markerfacecolor='#9C27B0',
        markeredgecolor='#4A148C',
        label=r'$\langle J \rangle$ (scanning rate)',
    )

    ax_left.set_xlabel('Numero de particulas $N$', fontsize=17)
    ax_left.set_ylabel(r'$J_{fresh}^{in}$ [particles/(m·s)]', color='#BF360C', fontsize=17)
    ax_right.set_ylabel(r'$\langle J \rangle$ [contactos/s]', color='#6A1B9A', fontsize=17)
    ax_left.set_xticks(N_values)
    ax_left.grid(True, alpha=0.3, linestyle='--')
    ax_left.tick_params(axis='both', labelsize=14)
    ax_right.tick_params(axis='both', labelsize=14)
    format_y_axis(ax_left, yerr=J_fresh_err)
    format_y_axis(ax_right, yerr=J_scan_err)

    handles = [left_line, right_line]
    labels = [h.get_label() for h in handles]
    ax_left.legend(handles, labels, loc='best', fontsize=13)

    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    return outpath


def compute_scanning_rate_summary(files_by_N):
    """Compute mean and SEM of scanning rate J for each N from cached entries."""
    N_values = sorted(files_by_N.keys())
    means = []
    errs = []
    for N in N_values:
        j_values = np.array([float(entry["cfc"]["J"]) for entry in files_by_N[N]], dtype=float)
        means.append(float(np.mean(j_values)))
        errs.append(float(np.std(j_values) / np.sqrt(len(j_values))) if len(j_values) > 0 else 0.0)
    return N_values, np.array(means, dtype=float), np.array(errs, dtype=float)


def plot_radial_profiles(entry=None, filepath=None, output_dir="graphics/output", dS=0.2):
    """
    Plot 3 radial profile curves: <rho_f^in>(S), |<vf_in>(S)|, J_in(S).
    """
    os.makedirs(output_dir, exist_ok=True)
    
    if entry is None:
        if filepath is None:
            print("  [1.4] No data provided.")
            return
        entry = load_analysis_file(filepath)
    
    metadata = entry["metadata"]
    N = int(metadata['N'])
    (
        S_centers,
        shell_areas,
        snapshot_count,
        count_sum,
        count_sq_sum,
        velocity_sample_count,
        vf_in_sum,
        vf_in_sq_sum,
    ), t_est = \
        extract_stationary_accumulators(entry, dS)
    rho_mean, _, v_mean, _, J_in = derive_profiles_from_stationary_samples(
        snapshot_count,
        shell_areas,
        count_sum,
        count_sq_sum,
        velocity_sample_count,
        vf_in_sum,
        vf_in_sq_sum,
    )

    outpaths = save_radial_profile_figures(
        N,
        S_centers,
        rho_mean,
        v_mean,
        J_in,
        output_dir,
    )
    for outpath in outpaths:
        print(f"  [1.4] Saved: {outpath}")
    print(f"  [1.4] N={N}: perfiles calculados usando snapshots con t>=t_est ({t_est}).")

    return S_centers, rho_mean, v_mean, J_in


def plot_radial_profiles_ensemble(entries, output_dir="graphics/output", dS=0.2):
    """Plot radial profiles averaged over all realizations for one N."""
    os.makedirs(output_dir, exist_ok=True)

    metadata = entries[0]["metadata"]
    N = int(metadata['N'])
    t_est = MANUAL_T_EST_BY_N.get(N)
    S_centers, rho_mean, rho_err, v_mean, v_err, J_mean, J_err = \
        aggregate_radial_profiles(entries, dS)

    outpaths = save_radial_profile_figures(
        N,
        S_centers,
        rho_mean,
        v_mean,
        J_mean,
        output_dir,
        rho_err=rho_err,
        v_err=v_err,
        J_err=J_err,
    )
    for outpath in outpaths:
        print(f"  [1.4] Saved: {outpath}")
    print(f"  [1.4] N={N}: promedio y barras SEM (std/sqrt(n)) calculados con t>=t_est ({t_est}).")

    return S_centers, rho_mean, v_mean, J_mean


def plot_at_S2_vs_N(files_by_N, output_dir="graphics/output", dS=0.2, s_range=NEAR_OBSTACLE_S_RANGE):
    """
    For the near-obstacle band S in [2,3], plot J_in, <rho_f^in>, and |<vf_in>| vs N.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    s_min, s_max = s_range

    N_values = sorted(files_by_N.keys())
    rho_at_S2 = []
    v_at_S2 = []
    J_at_S2 = []
    rho_err = []
    v_err = []
    j_err = []
    for N in N_values:
        stats = extract_profile_statistics_for_entries(files_by_N[N], dS)
        shell_mask = pick_shell_band_mask(stats["S_centers"], s_min, s_max)

        rho_mean_band, rho_err_band = reduce_shell_band(stats["rho_mean"], stats["rho_err"], shell_mask)
        v_mean_band, v_err_band = reduce_shell_band(np.abs(stats["v_mean"]), stats["v_err"], shell_mask)
        J_mean_band, J_err_band = reduce_shell_band(stats["J_mean"], stats["J_err"], shell_mask)

        rho_at_S2.append(rho_mean_band)
        v_at_S2.append(v_mean_band)
        J_at_S2.append(J_mean_band)
        rho_err.append(rho_err_band)
        v_err.append(v_err_band)
        j_err.append(J_err_band)

    outpaths = save_s2_vs_n_figures(
        N_values,
        rho_at_S2,
        v_at_S2,
        J_at_S2,
        output_dir,
        rho_err=rho_err,
        v_err=v_err,
        J_err=j_err,
    )
    for outpath in outpaths:
        print(f"  [1.4] Saved: {outpath}")

    return (
        np.array(N_values, dtype=float),
        np.array(rho_at_S2, dtype=float),
        np.array(v_at_S2, dtype=float),
        np.array(J_at_S2, dtype=float),
        np.array(rho_err, dtype=float),
        np.array(v_err, dtype=float),
        np.array(j_err, dtype=float),
    )


def plot_profiles_comparison_by_N(files_by_N, output_dir="graphics/output", dS=0.2, zoom_range=ZOOM_S_RANGE):
    """Plot radial profiles vs S for all N in one figure (+zoom)."""
    profile_by_N = collect_profiles_by_N(files_by_N, dS)
    full_path = save_profiles_by_n_figures(profile_by_N, output_dir, zoom_range=None)
    zoom_path = save_profiles_by_n_figures(profile_by_N, output_dir, zoom_range=zoom_range)
    print(f"  [1.4] Saved: {full_path}")
    print(f"  [1.4] Saved: {zoom_path}")


def plot_near_obstacle_vs_scanning_rate(files_by_N, output_dir="graphics/output", dS=0.2, s_range=NEAR_OBSTACLE_S_RANGE):
    """Compare near-obstacle J_fresh^in(N) with scanning-rate J(N)."""
    s2_data = plot_at_S2_vs_N(files_by_N, output_dir, dS=dS, s_range=s_range)
    N_values = s2_data[0].astype(int)
    J_fresh = s2_data[3]
    J_fresh_err = s2_data[6]

    N_scan, J_scan, J_scan_err = compute_scanning_rate_summary(files_by_N)
    scan_lookup = {N: (J_scan[idx], J_scan_err[idx]) for idx, N in enumerate(N_scan)}

    aligned_J_scan = []
    aligned_J_scan_err = []
    for N in N_values:
        mean_scan, err_scan = scan_lookup[int(N)]
        aligned_J_scan.append(mean_scan)
        aligned_J_scan_err.append(err_scan)

    compare_path = save_near_obstacle_vs_scanning_figure(
        N_values,
        J_fresh,
        J_fresh_err,
        np.array(aligned_J_scan, dtype=float),
        np.array(aligned_J_scan_err, dtype=float),
        output_dir,
    )
    print(f"  [1.4] Saved: {compare_path}")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python plot_radial_profiles.py <sim_file_or_directory>")
        sys.exit(1)
    
    path = sys.argv[1]
    
    if os.path.isfile(path):
        plot_radial_profiles(filepath=path)
    elif os.path.isdir(path):
        entries = load_analysis_entries(path)
        if not entries:
            print("No simulation files found.")
            sys.exit(1)
        files_by_N = group_entries_by_N(entries)
        selected_profile_ns = [N for N in sorted(files_by_N.keys()) if N in DISPLAY_PROFILE_N_VALUES]

        # Plot radial profiles for each N using all realizations.
        for N in selected_profile_ns:
            plot_radial_profiles_ensemble(files_by_N[N])

        plot_profiles_comparison_by_N(files_by_N)

        # Plot S∈[2,3] vs N and compare with scanning rate.
        if len(files_by_N) > 1:
            plot_near_obstacle_vs_scanning_rate(files_by_N)
    else:
        print(f"Error: '{path}' not found.")
        sys.exit(1)
