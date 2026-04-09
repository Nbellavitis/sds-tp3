"""
plot_scanning_rate.py - Inciso 1.2
Generates: Scanning rate <J>(N) with error bars.

The scanning rate J is obtained by linear interpolation of C_fc(t),
the cumulative count of fresh particles that contact the central obstacle.
J = slope of the linear fit of C_fc(t).

Usage:
    python graphics/plot_scanning_rate.py data/
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


def compute_cfc(snapshots):
    """
    Compute C_fc(t): cumulative count of fresh-to-center contacts.
    
    A "fresh contact with center" occurs when a particle transitions
    from FRESH to USED state between consecutive snapshots.
    
    Returns:
        times: array of snapshot times
        cfc: cumulative count of fresh->used transitions
    """
    times = []
    cfc_values = []
    cumulative = 0
    
    prev_states = None
    
    for t, particles in snapshots:
        current_states = {p['id']: p['state'] for p in particles}
        
        if prev_states is not None:
            # Count transitions F -> U (fresh particle hit the obstacle)
            for pid, state in current_states.items():
                if pid in prev_states and prev_states[pid] == 'F' and state == 'U':
                    cumulative += 1
        
        times.append(t)
        cfc_values.append(cumulative)
        prev_states = current_states
    
    return np.array(times), np.array(cfc_values)


def compute_scanning_rate(times, cfc):
    """
    Compute scanning rate J as the slope of linear fit of C_fc(t).
    Uses only the steady-state region (after transient).
    
    Returns: J (slope), intercept, r_squared
    """
    if len(times) < 2:
        return 0.0, 0.0, 0.0
    
    # Use the second half of the simulation for steady-state
    mid = len(times) // 4
    t_ss = times[mid:]
    c_ss = cfc[mid:]
    
    # Linear fit
    coeffs = np.polyfit(t_ss, c_ss, 1)
    J = coeffs[0]  # slope = scanning rate
    intercept = coeffs[1]
    
    # R-squared
    c_pred = np.polyval(coeffs, t_ss)
    ss_res = np.sum((c_ss - c_pred) ** 2)
    ss_tot = np.sum((c_ss - np.mean(c_ss)) ** 2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    
    return J, intercept, r_squared


def plot_scanning_rate(files_by_N=None, output_dir="graphics/output"):
    """
    Plot <J>(N) with error bars from multiple realizations.
    
    Args:
        files_by_N: dict {N: [(filepath, metadata, snapshots), ...]}
        output_dir: directory to save plots
    """
    os.makedirs(output_dir, exist_ok=True)
    
    if files_by_N is None:
        # Try loading from data directory
        data_dir = "data"
        files = sorted([os.path.join(data_dir, f) for f in os.listdir(data_dir)
                        if f.startswith('sim_') and f.endswith('.txt')])
        
        if not files:
            print("  [1.2] No simulation files found.")
            return
        
        files_by_N = {}
        for fpath in files:
            meta, snaps = parse_simulation_file(fpath)
            n = int(meta['N'])
            if n not in files_by_N:
                files_by_N[n] = []
            files_by_N[n].append((fpath, meta, snaps))
    
    N_values = sorted(files_by_N.keys())
    J_means = []
    J_stds = []
    
    for N in N_values:
        J_list = []
        for (fpath, meta, snaps) in files_by_N[N]:
            times, cfc = compute_cfc(snaps)
            J, _, r2 = compute_scanning_rate(times, cfc)
            J_list.append(J)
        
        J_means.append(np.mean(J_list))
        J_stds.append(np.std(J_list))
        print(f"  [1.2] N={N}: <J> = {np.mean(J_list):.4f} ± {np.std(J_list):.4f} contacts/s "
              f"({len(J_list)} runs)")
    
    # ── Plot <J>(N) ──────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 7))
    
    ax.errorbar(N_values, J_means, yerr=J_stds, fmt='s-', capsize=5,
                color='#9C27B0', ecolor='#CE93D8', markerfacecolor='#6A1B9A',
                markeredgecolor='#4A148C', markersize=8, linewidth=2,
                label='$\\langle J \\rangle \\pm \\sigma$')
    
    ax.set_xlabel('Número de partículas $N$', fontsize=14)
    ax.set_ylabel('Scanning rate $\\langle J \\rangle$ [contactos/s]', fontsize=14)
    ax.set_title('Inciso 1.2: Scanning rate $\\langle J \\rangle$ vs $N$',
                 fontsize=16, fontweight='bold')
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', labelsize=12)
    ax.ticklabel_format(axis='y', style='scientific', scilimits=(0, 3))
    
    plt.tight_layout()
    outpath = os.path.join(output_dir, "inciso_1_2_scanning_rate.png")
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [1.2] Saved: {outpath}")
    
    # ── Also plot C_fc(t) for the largest N as illustration ──────────────
    max_N = max(N_values)
    fig2, ax2 = plt.subplots(figsize=(10, 7))
    
    for fpath, meta, snaps in files_by_N[max_N]:
        times, cfc = compute_cfc(snaps)
        J, intercept, r2 = compute_scanning_rate(times, cfc)
        
        label_data = os.path.basename(fpath)
        ax2.plot(times, cfc, alpha=0.6, linewidth=1.5, label=f'{label_data}')
        
        # Plot linear fit
        t_fit = np.linspace(times[0], times[-1], 100)
        ax2.plot(t_fit, J * t_fit + intercept, '--', alpha=0.4, color='red')
    
    ax2.set_xlabel('Tiempo $t$ [s]', fontsize=14)
    ax2.set_ylabel('$C_{fc}(t)$ [contactos acumulados]', fontsize=14)
    ax2.set_title(f'Inciso 1.2: $C_{{fc}}(t)$ para $N={max_N}$\n'
                  f'(pendiente = scanning rate $J$)',
                  fontsize=16, fontweight='bold')
    ax2.legend(fontsize=10, loc='lower right')
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.tick_params(axis='both', labelsize=12)
    
    plt.tight_layout()
    outpath2 = os.path.join(output_dir, "inciso_1_2_cfc_curve.png")
    plt.savefig(outpath2, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [1.2] Saved: {outpath2}")


if __name__ == '__main__':
    plot_scanning_rate()
