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
- For S ≈ 2: plot J_in, <rho_f^in>, <vf_in> as functions of N

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


DISPLAY_PROFILE_N_VALUES = {100, 300, 500, 800}
STATIONARY_ACCUMULATORS_CACHE = {}


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
    """Stream one sim file and accumulate radial data only for snapshots with t >= t_est."""
    metadata = {}
    shell_edges = None
    S_centers = None
    shell_areas = None
    count_sum = None
    velocity_sum = None
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
                velocity_sum = np.zeros(len(S_centers), dtype=float)

            t_snapshot = float(trimmed.split()[1])
            include_snapshot = (t_est is None) or (t_snapshot >= t_est)

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
                count_sum[shell_idx] += 1.0
                velocity_sum[shell_idx] += vf_in

            if include_snapshot:
                snapshot_count += 1.0

    if S_centers is None:
        raise ValueError(f"No snapshots found in {filepath}")

    return S_centers, shell_areas, snapshot_count, count_sum, velocity_sum


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

    The mean profile is computed from pooled snapshot accumulators across all
    realizations, while the error bars still show per-realization dispersion.
    """
    stationary_data = [extract_stationary_accumulators(entry, dS) for entry in entries]
    accumulators = [item[0] for item in stationary_data]
    S_centers = accumulators[0][0]
    shell_areas = accumulators[0][1]

    pooled_snapshot_count = sum(acc[2] for acc in accumulators)
    pooled_count_sum = np.sum([acc[3] for acc in accumulators], axis=0)
    pooled_vf_in_sum = np.sum([acc[4] for acc in accumulators], axis=0)
    rho_mean, v_mean, J_mean = derive_profiles_from_accumulators(
        pooled_snapshot_count,
        shell_areas,
        pooled_count_sum,
        pooled_vf_in_sum,
    )

    per_run_profiles = [
        derive_profiles_from_accumulators(snapshot_count, shell_areas, count_sum, vf_in_sum)
        for _, _, snapshot_count, count_sum, vf_in_sum in accumulators
    ]
    rho_matrix = np.vstack([p[0] for p in per_run_profiles])
    v_matrix = np.vstack([p[1] for p in per_run_profiles])
    j_matrix = np.vstack([p[2] for p in per_run_profiles])

    return (
        S_centers,
        rho_mean,
        np.std(rho_matrix, axis=0, ddof=0),
        v_mean,
        np.std(v_matrix, axis=0, ddof=0),
        J_mean,
        np.std(j_matrix, axis=0, ddof=0),
    )


def derive_profiles_from_accumulators(snapshot_count, shell_areas, count_sum, vf_in_sum):
    """Build rho, v and J from pooled shell accumulators."""
    count_sum = np.array(count_sum, dtype=float)
    vf_in_sum = np.array(vf_in_sum, dtype=float)
    shell_areas = np.array(shell_areas, dtype=float)

    rho = np.zeros_like(count_sum, dtype=float)
    v = np.zeros_like(count_sum, dtype=float)

    if snapshot_count > 0:
        rho = count_sum / (float(snapshot_count) * shell_areas)

    nonzero = count_sum > 1e-12
    v[nonzero] = vf_in_sum[nonzero] / count_sum[nonzero]
    J = rho * np.abs(v)
    return rho, v, J


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

    ax.set_xlabel(xlabel, fontsize=13)
    ax.set_ylabel(ylabel, fontsize=13)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', labelsize=11)

    if xticks is not None:
        ax.set_xticks(xticks)

    if scientific_y:
        ax.ticklabel_format(axis='y', style='scientific', scilimits=(0, 3))

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
        r'$J_{in}$ [partículas/(m$^2$·s)]',
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
    """Save one file per observable for the shell nearest S≈2."""
    outputs = []

    rho_outpath = os.path.join(output_dir, "inciso_1_4_S2_rho_vs_N.png")
    save_profile_plot(
        N_values,
        rho_values,
        rho_outpath,
        'Número de partículas $N$',
        r'$\langle \rho_f^{in} \rangle$ [part/m²] en $S \approx 2$ m',
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

    v_outpath = os.path.join(output_dir, "inciso_1_4_S2_velocity_vs_N.png")
    save_profile_plot(
        N_values,
        v_values,
        v_outpath,
        'Número de partículas $N$',
        r'$|\langle v_f^{in} \rangle|$ [m/s] en $S \approx 2$ m',
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

    j_outpath = os.path.join(output_dir, "inciso_1_4_S2_flux_vs_N.png")
    save_profile_plot(
        N_values,
        J_values,
        j_outpath,
        'Número de partículas $N$',
        r'$J_{in}$ [part/(m²·s)] en $S \approx 2$ m',
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
    (S_centers, shell_areas, snapshot_count, count_sum, vf_in_sum), t_est = \
        extract_stationary_accumulators(entry, dS)
    rho_mean, v_mean, J_in = derive_profiles_from_accumulators(
        snapshot_count,
        shell_areas,
        count_sum,
        vf_in_sum,
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
    S_centers, rho_mean, rho_std, v_mean, v_std, J_mean, J_std = \
        aggregate_radial_profiles(entries, dS)

    outpaths = save_radial_profile_figures(
        N,
        S_centers,
        rho_mean,
        v_mean,
        J_mean,
        output_dir,
        rho_err=rho_std,
        v_err=v_std,
        J_err=J_std,
    )
    for outpath in outpaths:
        print(f"  [1.4] Saved: {outpath}")
    print(f"  [1.4] N={N}: promedio y desvío poblacional calculados con t>=t_est ({t_est}).")

    return S_centers, rho_mean, v_mean, J_mean


def plot_at_S2_vs_N(files_by_N, output_dir="graphics/output", dS=0.2):
    """
    For the shell at S ≈ 2m, plot J_in, <rho_f^in>, and |<vf_in>| as functions of N.
    
    Required by the enunciado: "Para la capa cercana a S=2, graficar Jin, 
    <ρ_f^in>, y <v_f^in> en función de N."
    """
    os.makedirs(output_dir, exist_ok=True)
    
    S_target = 2.0  # the enunciado specifies S ≈ 2
    
    N_values = sorted(files_by_N.keys())
    rho_at_S2 = []
    v_at_S2 = []
    J_at_S2 = []
    rho_err = []
    v_err = []
    j_err = []
    for N in N_values:
        per_run_profiles = []
        stationary_data = [extract_stationary_accumulators(entry, dS) for entry in files_by_N[N]]
        accumulators = [item[0] for item in stationary_data]
        S_centers = accumulators[0][0]
        shell_areas = accumulators[0][1]

        idx = np.where(S_centers >= S_target)[0][0]

        pooled_snapshot_count = sum(acc[2] for acc in accumulators)
        pooled_count_sum = np.sum([acc[3] for acc in accumulators], axis=0)
        pooled_vf_in_sum = np.sum([acc[4] for acc in accumulators], axis=0)
        rho_pooled, v_pooled, J_pooled = derive_profiles_from_accumulators(
            pooled_snapshot_count,
            shell_areas,
            pooled_count_sum,
            pooled_vf_in_sum,
        )

        for _, _, snapshot_count, count_sum, vf_in_sum in accumulators:
            per_run_profiles.append(
                derive_profiles_from_accumulators(
                    snapshot_count,
                    shell_areas,
                    count_sum,
                    vf_in_sum,
                )
            )

        rho_at_S2.append(rho_pooled[idx])
        v_at_S2.append(np.abs(v_pooled[idx]))
        J_at_S2.append(J_pooled[idx])
        rho_err.append(np.std([profile[0][idx] for profile in per_run_profiles], ddof=0))
        v_err.append(np.std([np.abs(profile[1][idx]) for profile in per_run_profiles], ddof=0))
        j_err.append(np.std([profile[2][idx] for profile in per_run_profiles], ddof=0))

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
        selected_files_by_N = {N: files_by_N[N] for N in selected_profile_ns}

        # Plot radial profiles for each N using all realizations.
        for N in selected_profile_ns:
            plot_radial_profiles_ensemble(files_by_N[N])
        
        # Plot S≈2 vs N
        if len(selected_files_by_N) > 1:
            plot_at_S2_vs_N(selected_files_by_N)
    else:
        print(f"Error: '{path}' not found.")
        sys.exit(1)
