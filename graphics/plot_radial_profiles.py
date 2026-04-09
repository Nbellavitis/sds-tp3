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


def compute_radial_profiles(snapshots, metadata, dS=0.2):
    """
    Compute radial profiles for fresh inward-moving particles.
    
    S = radial distance from center (origin).
    Shells from S=0 to S=R_enclosure.
    Only particles within valid center range (r0+r to R_enc-r) contribute.
    
    For each snapshot (steady state only), for each shell:
    - Count fresh particles with R·v < 0
    - Compute density = count / shell_area
    - Compute mean radial velocity: vf_in = (R·v)/|R| (negative = inward)
    
    Returns:
        S_centers: array of shell center distances
        rho_mean, rho_std: density <rho_f^in>(S) ± error
        v_mean, v_std: radial velocity <vf_in>(S) ± error (not normalized)
        J_in: flux = rho * |v|
    """
    R_enclosure = float(metadata['R_enclosure'])
    r0 = float(metadata['r0'])
    r_particle = float(metadata['r'])
    
    # S ranges from 0 to R_enclosure
    # But particle centers can only be at S in [r0+r, R_enc-r]
    S_min = 0.0
    S_max = R_enclosure
    
    shell_edges = np.arange(S_min, S_max + dS, dS)
    n_shells = len(shell_edges) - 1
    S_centers = (shell_edges[:-1] + shell_edges[1:]) / 2.0
    
    # Shell areas
    shell_areas = np.array([
        np.pi * (shell_edges[k+1]**2 - shell_edges[k]**2)
        for k in range(n_shells)
    ])
    
    # Only use steady-state snapshots (last 50%)
    ss_start = len(snapshots) // 2
    ss_snapshots = snapshots[ss_start:]
    n_frames = len(ss_snapshots)
    
    if n_frames == 0:
        return S_centers, np.zeros(n_shells), np.zeros(n_shells), \
               np.zeros(n_shells), np.zeros(n_shells), np.zeros(n_shells)
    
    # Per-frame accumulation
    count_per_frame = np.zeros((n_frames, n_shells))
    v_sum_per_frame = np.zeros((n_frames, n_shells))
    
    for frame_idx, (t, particles) in enumerate(ss_snapshots):
        for p in particles:
            # Filter: only FRESH particles
            if p['state'] != 'F':
                continue
            
            x, y = p['x'], p['y']
            vx, vy = p['vx'], p['vy']
            
            S = np.sqrt(x * x + y * y)  # distance from center
            
            # Filter: R · v < 0 (velocity pointing toward center)
            rdotv = x * vx + y * vy
            if rdotv >= 0:
                continue
            
            # Shell index
            shell_idx = int((S - S_min) / dS)
            if shell_idx < 0 or shell_idx >= n_shells:
                continue
            
            # Radial velocity: vf_in = (R · v) / |R|
            # This is NEGATIVE when inward (which we want)
            vf_in = rdotv / S if S > 1e-10 else 0.0
            
            count_per_frame[frame_idx, shell_idx] += 1
            v_sum_per_frame[frame_idx, shell_idx] += vf_in
    
    # Density per frame: count / area
    rho_per_frame = np.zeros((n_frames, n_shells))
    v_per_frame = np.zeros((n_frames, n_shells))
    
    for k in range(n_shells):
        rho_per_frame[:, k] = count_per_frame[:, k] / shell_areas[k]
        for f in range(n_frames):
            if count_per_frame[f, k] > 0:
                v_per_frame[f, k] = v_sum_per_frame[f, k] / count_per_frame[f, k]
    
    # Average over frames
    rho_mean = np.mean(rho_per_frame, axis=0)
    rho_std = np.std(rho_per_frame, axis=0) / np.sqrt(n_frames)
    
    v_mean = np.mean(v_per_frame, axis=0)  # negative (inward)
    v_std = np.std(v_per_frame, axis=0) / np.sqrt(n_frames)
    
    # J_in = <rho> * |<v>|
    J_in = rho_mean * np.abs(v_mean)
    
    return S_centers, rho_mean, rho_std, v_mean, v_std, J_in


def plot_radial_profiles(metadata=None, snapshots=None, filepath=None,
                         output_dir="graphics/output", dS=0.2):
    """
    Plot 3 radial profile curves: <rho_f^in>(S), |<vf_in>(S)|, J_in(S).
    """
    os.makedirs(output_dir, exist_ok=True)
    
    if metadata is None or snapshots is None:
        if filepath is None:
            print("  [1.4] No data provided.")
            return
        metadata, snapshots, _ = parse_simulation_file(filepath)
    
    N = int(metadata['N'])
    r0 = float(metadata['r0'])
    r_particle = float(metadata['r'])
    
    S_centers, rho_mean, rho_std, v_mean, v_std, J_in = \
        compute_radial_profiles(snapshots, metadata, dS)
    
    # Mask: only plot where there's valid physical range
    valid = S_centers >= (r0 + r_particle - dS)
    
    # ── Three-panel figure ────────────────────────────────────────────────
    fig, axes = plt.subplots(3, 1, figsize=(10, 18), sharex=True)
    
    # --- Panel 1: Density <rho_f^in>(S) ---
    ax1 = axes[0]
    ax1.errorbar(S_centers[valid], rho_mean[valid], yerr=rho_std[valid],
                 fmt='o-', capsize=3, color='#1B5E20', markerfacecolor='#4CAF50',
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
    ax2.errorbar(S_centers[valid], np.abs(v_mean[valid]), yerr=v_std[valid],
                 fmt='s-', capsize=3, color='#0D47A1', markerfacecolor='#2196F3',
                 markeredgecolor='#0D47A1', markersize=3, linewidth=1.2,
                 label=r'$|\langle v_f^{in} \rangle(S)|$')
    ax2.set_ylabel(r'$|\langle v_f^{in} \rangle|$ [m/s]', fontsize=13)
    ax2.legend(fontsize=12)
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.tick_params(axis='both', labelsize=11)
    
    # --- Panel 3: Flux J_in(S) = <rho> * |<v>| ---
    ax3 = axes[2]
    ax3.plot(S_centers[valid], J_in[valid], 'D-', color='#BF360C',
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
    
    for N in N_values:
        rho_list = []
        v_list = []
        j_list = []
        
        for entry in files_by_N[N]:
            fpath, meta, snaps = entry[0], entry[1], entry[2]
            S_centers, rho, _, v, _, J = compute_radial_profiles(snaps, meta, dS)
            
            # Find shell closest to S=2
            idx = np.argmin(np.abs(S_centers - S_target))
            rho_list.append(rho[idx])
            v_list.append(np.abs(v[idx]))
            j_list.append(J[idx])
        
        rho_at_S2.append(np.mean(rho_list))
        v_at_S2.append(np.mean(v_list))
        J_at_S2.append(np.mean(j_list))
        rho_err.append(np.std(rho_list))
        v_err.append(np.std(v_list))
    
    # ── Plot ──────────────────────────────────────────────────────────────
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 15), sharex=True)
    
    ax1.errorbar(N_values, rho_at_S2, yerr=rho_err, fmt='o-', capsize=5,
                 color='#1B5E20', markerfacecolor='#4CAF50', markersize=8,
                 linewidth=2, label=r'$\langle \rho_f^{in} \rangle$ at $S \approx 2$ m')
    ax1.set_ylabel(r'$\langle \rho_f^{in} \rangle$ [part/m²]', fontsize=13)
    ax1.set_title(f'Inciso 1.4: Perfiles radiales en $S \\approx {S_target}$ m vs $N$',
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
    
    ax3.plot(N_values, J_at_S2, 'D-', color='#BF360C',
             markerfacecolor='#FF5722', markersize=8, linewidth=2,
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
        metadata, snapshots, _ = parse_simulation_file(path)
        plot_radial_profiles(metadata, snapshots, path)
    elif os.path.isdir(path):
        files = sorted(glob.glob(os.path.join(path, 'sim_*.txt')))
        if not files:
            print("No simulation files found.")
            sys.exit(1)
        
        files_by_N = {}
        for fpath in files:
            meta, snaps, events = parse_simulation_file(fpath)
            n = int(meta['N'])
            if n not in files_by_N:
                files_by_N[n] = []
            files_by_N[n].append((fpath, meta, snaps, events))
        
        # Plot radial profiles for each N (first realization)
        for N in sorted(files_by_N.keys()):
            fpath, meta, snaps, _ = files_by_N[N][0]
            plot_radial_profiles(meta, snaps, fpath)
        
        # Plot S≈2 vs N
        if len(files_by_N) > 1:
            plot_at_S2_vs_N(files_by_N)
    else:
        print(f"Error: '{path}' not found.")
        sys.exit(1)
