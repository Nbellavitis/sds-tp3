"""
plot_fraction_used.py - Inciso 1.3
Generates: Temporal evolution of the fraction of used particles F_u(t).

Shows how the fraction of "used" (violet) particles evolves over time
and identifies the steady-state regime.

Usage:
    python graphics/plot_fraction_used.py data/sim_100N_20260410_1530.txt
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


def compute_fraction_used(snapshots, N):
    """
    Compute F_u(t) = fraction of particles in USED state at each snapshot.
    
    Returns:
        times: array of snapshot times
        fu: array of fraction used at each time
    """
    times = []
    fu = []
    
    for t, particles in snapshots:
        n_used = sum(1 for p in particles if p['state'] == 'U')
        times.append(t)
        fu.append(n_used / N)
    
    return np.array(times), np.array(fu)


def plot_fraction_used(metadata=None, snapshots=None, filepath=None,
                       output_dir="graphics/output"):
    """
    Plot F_u(t): fraction of used particles over time.
    Shows transient and steady-state regimes.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    if metadata is None or snapshots is None:
        if filepath is None:
            print("  [1.3] No data provided.")
            return
        metadata, snapshots = parse_simulation_file(filepath)
    
    N = int(metadata['N'])
    times, fu = compute_fraction_used(snapshots, N)
    
    if len(times) == 0:
        print("  [1.3] No snapshots found.")
        return
    
    # ── Detect steady state ───────────────────────────────────────────────
    # Moving average for smoothing
    window = max(1, len(fu) // 20)
    if window > 1:
        fu_smooth = np.convolve(fu, np.ones(window) / window, mode='valid')
        t_smooth = times[:len(fu_smooth)]
    else:
        fu_smooth = fu
        t_smooth = times
    
    # Estimate steady-state value (mean of last 25% of data)
    ss_start = int(0.75 * len(fu))
    fu_ss = np.mean(fu[ss_start:])
    fu_ss_std = np.std(fu[ss_start:])
    
    # ── Plot ──────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 7))
    
    ax.plot(times, fu, alpha=0.4, color='#7B1FA2', linewidth=0.8,
            label='$F_u(t)$ instantáneo')
    
    if window > 1:
        ax.plot(t_smooth, fu_smooth, color='#4A148C', linewidth=2.0,
                label=f'$F_u(t)$ promedio móvil (ventana={window})')
    
    # Steady-state band
    ax.axhline(y=fu_ss, color='#E91E63', linestyle='--', linewidth=1.5,
               label=f'Estacionario: $F_u \\approx {fu_ss:.3f} \\pm {fu_ss_std:.3f}$')
    ax.axhspan(fu_ss - fu_ss_std, fu_ss + fu_ss_std, alpha=0.1, color='#E91E63')
    
    ax.set_xlabel('Tiempo $t$ [s]', fontsize=14)
    ax.set_ylabel('Fracción de partículas usadas $F_u(t)$', fontsize=14)
    ax.set_title(f'Inciso 1.3: Evolución temporal de $F_u(t)$ ($N={N}$)',
                 fontsize=16, fontweight='bold')
    ax.legend(fontsize=12, loc='lower right')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', labelsize=12)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlim(times[0], times[-1])
    
    plt.tight_layout()
    outpath = os.path.join(output_dir, f"inciso_1_3_fraction_used_N{N}.png")
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [1.3] Saved: {outpath}")
    print(f"  [1.3] Steady-state F_u = {fu_ss:.4f} ± {fu_ss_std:.4f}")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python plot_fraction_used.py <sim_file>")
        sys.exit(1)
    
    filepath = sys.argv[1]
    metadata, snapshots = parse_simulation_file(filepath)
    plot_fraction_used(metadata, snapshots, filepath)
