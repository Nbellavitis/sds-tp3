"""
plot_radial_profiles.py - Inciso 1.4
Generates: Radial profiles of density, velocity, and flux for "fresh inward" particles.

Divides the annular space between obstacle and outer wall into concentric shells
of width dS = 0.2 m. Filters only FRESH particles with velocity pointing toward
the center (R_j · v_j < 0). Computes:
  - <rho_f^in>(S): density of fresh inward particles per shell
  - <v_f^in>(S): normalized velocity of fresh inward particles per shell
  - J_in(S) = <rho_f^in> * <v_f^in>: inward flux

Usage:
    python graphics/plot_radial_profiles.py data/sim_100N_20260410_1530.txt
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def parse_simulation_file(filepath):
    """Parse simulation file and return metadata + snapshots."""
    metadata = {}
    snapshots = []
    
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
                particles.append({
                    'id': int(parts[0]),
                    'x': float(parts[1]),
                    'y': float(parts[2]),
                    'vx': float(parts[3]),
                    'vy': float(parts[4]),
                    'state': parts[5]
                })
            snapshots.append((time, particles))
        i += 1
    
    return metadata, snapshots


def compute_radial_profiles(snapshots, metadata, dS=0.2):
    """
    Compute radial profiles for fresh inward-moving particles.
    
    Only considers particles in steady state (last 50% of simulation).
    
    Args:
        snapshots: list of (time, particles_data)
        metadata: dict with simulation parameters
        dS: shell width [m]
    
    Returns:
        r_centers: array of shell center radii
        rho_mean, rho_std: density profile with error
        v_mean, v_std: velocity profile with error
        J_in: flux profile (rho * v)
    """
    R_enclosure = float(metadata['R_enclosure'])
    r0 = float(metadata['r0'])
    r_particle = float(metadata['r'])
    
    # Radial range: from obstacle surface + particle radius to enclosure - particle radius
    r_min = r0 + r_particle
    r_max = R_enclosure - r_particle
    
    # Shell edges
    shell_edges = np.arange(r_min, r_max + dS, dS)
    n_shells = len(shell_edges) - 1
    r_centers = (shell_edges[:-1] + shell_edges[1:]) / 2.0
    
    # Only use steady-state snapshots (last 50%)
    ss_start = len(snapshots) // 2
    ss_snapshots = snapshots[ss_start:]
    
    n_frames = len(ss_snapshots)
    
    # Accumulate per-frame profiles
    rho_profiles = np.zeros((n_frames, n_shells))
    v_profiles = np.zeros((n_frames, n_shells))
    count_profiles = np.zeros((n_frames, n_shells))
    
    for frame_idx, (t, particles) in enumerate(ss_snapshots):
        for p in particles:
            # Filter: only FRESH particles
            if p['state'] != 'F':
                continue
            
            x, y = p['x'], p['y']
            vx, vy = p['vx'], p['vy']
            
            r = np.sqrt(x * x + y * y)
            
            # Filter: velocity pointing toward center (R · v < 0)
            rdotv = x * vx + y * vy
            if rdotv >= 0:
                continue
            
            # Find which shell this particle is in
            shell_idx = int((r - r_min) / dS)
            if shell_idx < 0 or shell_idx >= n_shells:
                continue
            
            # Compute radial inward velocity component
            v_radial_inward = -rdotv / r  # positive when moving inward
            v_magnitude = np.sqrt(vx * vx + vy * vy)
            
            count_profiles[frame_idx, shell_idx] += 1
            v_profiles[frame_idx, shell_idx] += v_radial_inward
    
    # Compute density: number / shell area
    for k in range(n_shells):
        r_inner = shell_edges[k]
        r_outer = shell_edges[k + 1]
        shell_area = np.pi * (r_outer ** 2 - r_inner ** 2)
        rho_profiles[:, k] = count_profiles[:, k] / shell_area
    
    # Compute average velocity per shell per frame
    for frame_idx in range(n_frames):
        for k in range(n_shells):
            if count_profiles[frame_idx, k] > 0:
                v_profiles[frame_idx, k] /= count_profiles[frame_idx, k]
    
    # Average across frames
    rho_mean = np.mean(rho_profiles, axis=0)
    rho_std = np.std(rho_profiles, axis=0) / np.sqrt(n_frames)  # SEM
    
    v_mean = np.mean(v_profiles, axis=0)
    v_std = np.std(v_profiles, axis=0) / np.sqrt(n_frames)
    
    # Normalize velocity by v0
    v0 = float(metadata.get('v0', 1.0))
    v_mean_norm = v_mean / v0
    v_std_norm = v_std / v0
    
    # Inward flux
    J_in = rho_mean * v_mean
    
    return r_centers, rho_mean, rho_std, v_mean_norm, v_std_norm, J_in


def plot_radial_profiles(metadata=None, snapshots=None, filepath=None,
                         output_dir="graphics/output", dS=0.2):
    """
    Plot radial profiles: density, velocity, and flux.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    if metadata is None or snapshots is None:
        if filepath is None:
            print("  [1.4] No data provided.")
            return
        metadata, snapshots = parse_simulation_file(filepath)
    
    N = int(metadata['N'])
    
    r_centers, rho_mean, rho_std, v_mean, v_std, J_in = \
        compute_radial_profiles(snapshots, metadata, dS)
    
    # ── Three-panel figure ────────────────────────────────────────────────
    fig, axes = plt.subplots(3, 1, figsize=(10, 18), sharex=True)
    
    # --- Panel 1: Density ---
    ax1 = axes[0]
    ax1.errorbar(r_centers, rho_mean, yerr=rho_std, fmt='o-', capsize=3,
                 color='#1B5E20', markerfacecolor='#4CAF50',
                 markeredgecolor='#1B5E20', markersize=4, linewidth=1.5,
                 label=r'$\langle \rho_f^{in} \rangle(S)$')
    ax1.set_ylabel(r'Densidad $\langle \rho_f^{in} \rangle$ [partículas/m²]',
                   fontsize=13)
    ax1.set_title(f'Inciso 1.4: Perfiles radiales de partículas frescas entrantes ($N={N}$, $dS={dS}$ m)',
                  fontsize=15, fontweight='bold')
    ax1.legend(fontsize=12)
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.tick_params(axis='both', labelsize=11)
    ax1.ticklabel_format(axis='y', style='scientific', scilimits=(0, 3))
    
    # --- Panel 2: Normalized velocity ---
    ax2 = axes[1]
    ax2.errorbar(r_centers, v_mean, yerr=v_std, fmt='s-', capsize=3,
                 color='#0D47A1', markerfacecolor='#2196F3',
                 markeredgecolor='#0D47A1', markersize=4, linewidth=1.5,
                 label=r'$\langle v_f^{in} \rangle / v_0$')
    ax2.set_ylabel(r'Velocidad normalizada $\langle v_f^{in} \rangle / v_0$',
                   fontsize=13)
    ax2.legend(fontsize=12)
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.tick_params(axis='both', labelsize=11)
    
    # --- Panel 3: Flux ---
    ax3 = axes[2]
    ax3.plot(r_centers, J_in, 'D-', color='#BF360C', markerfacecolor='#FF5722',
             markeredgecolor='#BF360C', markersize=4, linewidth=1.5,
             label=r'$J_{in}(S) = \langle \rho_f^{in} \rangle \cdot \langle v_f^{in} \rangle$')
    ax3.set_xlabel('Distancia radial $S$ [m]', fontsize=13)
    ax3.set_ylabel(r'Flujo entrante $J_{in}$ [partículas/(m·s)]', fontsize=13)
    ax3.legend(fontsize=12)
    ax3.grid(True, alpha=0.3, linestyle='--')
    ax3.tick_params(axis='both', labelsize=11)
    ax3.ticklabel_format(axis='y', style='scientific', scilimits=(0, 3))
    
    plt.tight_layout()
    outpath = os.path.join(output_dir, f"inciso_1_4_radial_profiles_N{N}.png")
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [1.4] Saved: {outpath}")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python plot_radial_profiles.py <sim_file>")
        sys.exit(1)
    
    filepath = sys.argv[1]
    metadata, snapshots = parse_simulation_file(filepath)
    plot_radial_profiles(metadata, snapshots, filepath)
