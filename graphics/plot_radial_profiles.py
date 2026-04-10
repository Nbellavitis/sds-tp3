"""
plot_radial_profiles.py - Inciso 1.4
Generates: Radial profiles for fresh inward-moving particles.

Per the enunciado:
- S = distance from center of fixed particle (= origin)
- dS = 0.2 m concentric shell width
- Select only fresh particles Pf^in with Rj · vj < 0 (velocity toward center)
- Compute:
  - <rho_f^in>(S): density = count / shell_area, averaged over time & realizations
  - <vf_in>(S): radial velocity vf_in_j = (Rj · vj) / |Rj|, averaged
  - J_in(S) = <rho_f^in>(S) * |<vf_in>(S)|
- Plot 3 curves: <rho_f^in>(S), |<vf_in>(S)|, J_in(S)
- For S ≈ 2: plot J_in, <rho_f^in>, <vf_in> as functions of N

Usage:
    python graphics/plot_radial_profiles.py data/sim_300N_file.txt
    python graphics/plot_radial_profiles.py data/
"""

import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from analysis_cache import group_entries_by_N, load_analysis_entries, load_analysis_file


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
    
    S = radial distance from center (origin).
    Shells cover the physically accessible range of particle centers:
    S in [r0 + r, R_enclosure - r].

    For each recorded event interval, for each shell:
    - Count fresh particles with R·v < 0
    - Compute density = count / shell_area
    - Compute mean radial velocity: vf_in = (R·v)/|R| (negative = inward)

    The averages are time-weighted using the interval between consecutive
    snapshots. This avoids biasing the profiles toward periods with many events.
    
    Returns:
        S_centers: array of shell center distances
        rho_mean: density <rho_f^in>(S)
        v_mean: radial velocity <vf_in>(S)
        J_in: flux = rho * |v|
    """
    shell_edges, S_centers, shell_areas = build_shell_geometry(metadata, dS)
    n_shells = len(S_centers)

    if len(snapshots) < 2:
        return S_centers, np.zeros(n_shells), np.zeros(n_shells)

    times = np.array([t for t, _ in snapshots], dtype=float)
    durations = np.diff(times)
    total_time = np.sum(durations)

    if total_time <= 0:
        return S_centers, np.zeros(n_shells), np.zeros(n_shells)

    count_time = np.zeros(n_shells)
    velocity_time_sum = np.zeros(n_shells)
    particle_time = np.zeros(n_shells)

    for (t, particles), dt in zip(snapshots[:-1], durations):
        if dt <= 0:
            continue

        counts = np.zeros(n_shells)
        velocity_sum = np.zeros(n_shells)

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
            counts[shell_idx] += 1
            velocity_sum[shell_idx] += vf_in

        count_time += dt * counts
        velocity_time_sum += dt * velocity_sum
        particle_time += dt * counts

    rho_mean = count_time / (total_time * shell_areas)
    v_mean = np.zeros(n_shells)
    nonzero = particle_time > 1e-12
    v_mean[nonzero] = velocity_time_sum[nonzero] / particle_time[nonzero]
    J_in = rho_mean * np.abs(v_mean)

    return S_centers, rho_mean, v_mean, J_in


def aggregate_radial_profiles(entries, dS=0.2):
    """Average radial profiles over multiple realizations with per-run dispersion."""
    profiles = [extract_radial_profiles(entry) for entry in entries]

    S_centers = profiles[0][0]
    rho_matrix = np.vstack([p[1] for p in profiles])
    v_matrix = np.vstack([p[2] for p in profiles])
    j_matrix = np.vstack([p[3] for p in profiles])

    return (
        S_centers,
        np.mean(rho_matrix, axis=0),
        np.std(rho_matrix, axis=0),
        np.mean(v_matrix, axis=0),
        np.std(v_matrix, axis=0),
        np.mean(j_matrix, axis=0),
        np.std(j_matrix, axis=0),
    )


def extract_radial_profiles(entry):
    """Return cached radial profile arrays for one realization."""
    radial = entry["radial_profiles"]
    return (
        np.array(radial["S_centers"], dtype=float),
        np.array(radial["rho"], dtype=float),
        np.array(radial["v"], dtype=float),
        np.array(radial["J"], dtype=float),
    )


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
    S_centers, rho_mean, v_mean, J_in = extract_radial_profiles(entry)
    
    # ── Three-panel figure ────────────────────────────────────────────────
    fig, axes = plt.subplots(3, 1, figsize=(10, 18), sharex=True)
    
    # --- Panel 1: Density <rho_f^in>(S) ---
    ax1 = axes[0]
    ax1.plot(S_centers, rho_mean, 'o-', color='#1B5E20', markerfacecolor='#4CAF50',
             markeredgecolor='#1B5E20', markersize=3, linewidth=1.2,
             label=r'$\langle \rho_f^{in} \rangle(S)$')
    ax1.set_ylabel(r'$\langle \rho_f^{in} \rangle$ [partículas/m²]', fontsize=13)
    ax1.set_title(f'Inciso 1.4: Perfiles radiales — partículas frescas entrantes '
                  f'($N={N}$, $dS={dS}$ m)', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=12)
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.tick_params(axis='both', labelsize=11)
    ax1.ticklabel_format(axis='y', style='scientific', scilimits=(0, 3))
    
    # --- Panel 2: |<vf_in>(S)| ---
    ax2 = axes[1]
    ax2.plot(S_centers, np.abs(v_mean), 's-', color='#0D47A1', markerfacecolor='#2196F3',
             markeredgecolor='#0D47A1', markersize=3, linewidth=1.2,
             label=r'$|\langle v_f^{in} \rangle(S)|$')
    ax2.set_ylabel(r'$|\langle v_f^{in} \rangle|$ [m/s]', fontsize=13)
    ax2.legend(fontsize=12)
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.tick_params(axis='both', labelsize=11)
    
    # --- Panel 3: Flux J_in(S) = <rho> * |<v>| ---
    ax3 = axes[2]
    ax3.plot(S_centers, J_in, 'D-', color='#BF360C',
             markerfacecolor='#FF5722', markeredgecolor='#BF360C',
             markersize=3, linewidth=1.2,
             label=r'$J_{in}(S) = \langle \rho_f^{in} \rangle \cdot |\langle v_f^{in} \rangle|$')
    ax3.set_xlabel('Distancia radial $S$ [m]', fontsize=13)
    ax3.set_ylabel(r'$J_{in}$ [partículas/(m$^2$·s)]', fontsize=13)
    ax3.legend(fontsize=12)
    ax3.grid(True, alpha=0.3, linestyle='--')
    ax3.tick_params(axis='both', labelsize=11)
    ax3.ticklabel_format(axis='y', style='scientific', scilimits=(0, 3))
    
    plt.tight_layout()
    outpath = os.path.join(output_dir, f"inciso_1_4_radial_profiles_N{N}.png")
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [1.4] Saved: {outpath}")
    
    return S_centers, rho_mean, v_mean, J_in


def plot_radial_profiles_ensemble(entries, output_dir="graphics/output", dS=0.2):
    """Plot radial profiles averaged over all realizations for one N."""
    os.makedirs(output_dir, exist_ok=True)

    metadata = entries[0]["metadata"]
    N = int(metadata['N'])
    n_runs = len(entries)

    S_centers, rho_mean, rho_std, v_mean, v_std, J_mean, J_std = \
        aggregate_radial_profiles(entries, dS)

    fig, axes = plt.subplots(3, 1, figsize=(10, 18), sharex=True)

    ax1 = axes[0]
    ax1.errorbar(S_centers, rho_mean, yerr=rho_std, fmt='o-', capsize=3,
                 color='#1B5E20', markerfacecolor='#4CAF50',
                 markeredgecolor='#1B5E20', markersize=3, linewidth=1.2,
                 label=rf'$\langle \rho_f^{{in}} \rangle(S)$ ({n_runs} realizaciones)')
    ax1.set_ylabel(r'$\langle \rho_f^{in} \rangle$ [partículas/m²]', fontsize=13)
    ax1.set_title(f'Inciso 1.4: Perfiles radiales — partículas frescas entrantes '
                  f'($N={N}$, $dS={dS}$ m)\npromedio sobre {n_runs} realizaciones',
                  fontsize=14, fontweight='bold')
    ax1.legend(fontsize=12)
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.tick_params(axis='both', labelsize=11)
    ax1.ticklabel_format(axis='y', style='scientific', scilimits=(0, 3))

    ax2 = axes[1]
    ax2.errorbar(S_centers, np.abs(v_mean), yerr=v_std, fmt='s-', capsize=3,
                 color='#0D47A1', markerfacecolor='#2196F3',
                 markeredgecolor='#0D47A1', markersize=3, linewidth=1.2,
                 label=r'$|\langle v_f^{in} \rangle(S)|$')
    ax2.set_ylabel(r'$|\langle v_f^{in} \rangle|$ [m/s]', fontsize=13)
    ax2.legend(fontsize=12)
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.tick_params(axis='both', labelsize=11)

    ax3 = axes[2]
    ax3.errorbar(S_centers, J_mean, yerr=J_std, fmt='D-', capsize=3,
                 color='#BF360C', markerfacecolor='#FF5722',
                 markeredgecolor='#BF360C', markersize=3, linewidth=1.2,
                 label=r'$J_{in}(S)$')
    ax3.set_xlabel('Distancia radial $S$ [m]', fontsize=13)
    ax3.set_ylabel(r'$J_{in}$ [partículas/(m$^2$·s)]', fontsize=13)
    ax3.legend(fontsize=12)
    ax3.grid(True, alpha=0.3, linestyle='--')
    ax3.tick_params(axis='both', labelsize=11)
    ax3.ticklabel_format(axis='y', style='scientific', scilimits=(0, 3))

    plt.tight_layout()
    outpath = os.path.join(output_dir, f"inciso_1_4_radial_profiles_N{N}.png")
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [1.4] Saved: {outpath}")

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
    selected_center = None
    
    for N in N_values:
        rho_list = []
        v_list = []
        j_list = []
        
        for entry in files_by_N[N]:
            S_centers, rho, v, J = extract_radial_profiles(entry)

            # Pick the first physically valid shell, which is the one nearest S=2.
            idx = np.where(S_centers >= S_target)[0][0]
            selected_center = S_centers[idx]
            rho_list.append(rho[idx])
            v_list.append(np.abs(v[idx]))
            j_list.append(J[idx])
        
        rho_at_S2.append(np.mean(rho_list))
        v_at_S2.append(np.mean(v_list))
        J_at_S2.append(np.mean(j_list))
        rho_err.append(np.std(rho_list))
        v_err.append(np.std(v_list))
        j_err.append(np.std(j_list))
    
    # ── Plot ──────────────────────────────────────────────────────────────
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 15), sharex=True)
    
    ax1.errorbar(N_values, rho_at_S2, yerr=rho_err, fmt='o-', capsize=5,
                 color='#1B5E20', markerfacecolor='#4CAF50', markersize=8,
                 linewidth=2, label=r'$\langle \rho_f^{in} \rangle$ at $S \approx 2$ m')
    ax1.set_ylabel(r'$\langle \rho_f^{in} \rangle$ [part/m²]', fontsize=13)
    ax1.set_title(f'Inciso 1.4: Perfiles radiales en la capa mas cercana a '
                  f'$S={S_target}$ m ($S_c={selected_center:.1f}$ m) vs $N$',
                  fontsize=15, fontweight='bold')
    ax1.legend(fontsize=12)
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.ticklabel_format(axis='y', style='scientific', scilimits=(0, 3))
    
    ax2.errorbar(N_values, v_at_S2, yerr=v_err, fmt='s-', capsize=5,
                 color='#0D47A1', markerfacecolor='#2196F3', markersize=8,
                 linewidth=2, label=r'$|\langle v_f^{in} \rangle|$ at $S \approx 2$ m')
    ax2.set_ylabel(r'$|\langle v_f^{in} \rangle|$ [m/s]', fontsize=13)
    ax2.legend(fontsize=12)
    ax2.grid(True, alpha=0.3, linestyle='--')
    
    ax3.errorbar(N_values, J_at_S2, yerr=j_err, fmt='D-', capsize=5,
             color='#BF360C', markerfacecolor='#FF5722', markersize=8, linewidth=2,
             label=r'$J_{in}$ at $S \approx 2$ m')
    ax3.set_xlabel('Número de partículas $N$', fontsize=13)
    ax3.set_ylabel(r'$J_{in}$ [part/(m²·s)]', fontsize=13)
    ax3.legend(fontsize=12)
    ax3.grid(True, alpha=0.3, linestyle='--')
    ax3.ticklabel_format(axis='y', style='scientific', scilimits=(0, 3))
    
    plt.tight_layout()
    outpath = os.path.join(output_dir, "inciso_1_4_S2_vs_N.png")
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
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
        
        # Plot radial profiles for each N using all realizations.
        for N in sorted(files_by_N.keys()):
            plot_radial_profiles_ensemble(files_by_N[N])
        
        # Plot S≈2 vs N
        if len(files_by_N) > 1:
            plot_at_S2_vs_N(files_by_N)
    else:
        print(f"Error: '{path}' not found.")
        sys.exit(1)
