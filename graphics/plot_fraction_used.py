"""
plot_fraction_used.py - Inciso 1.3
Generates: Temporal evolution of the fraction of used particles F_u(t).

Per the enunciado:
- "Estudiar la evolución temporal de la fracción de partículas usadas: Fu(t) = Nu(t)/N"
- "Reportar tiempo al estacionario y valor del estacionario alcanzado (Fest) en función de N"

This script generates:
1. F_u(t) vs t plot for each N (shows transient and steady-state)
2. F_est vs N plot (steady-state fraction as a function of N)
3. t_est vs N plot (time to reach steady state)

Usage:
    python graphics/plot_fraction_used.py data/sim_300N_*.txt
    python graphics/plot_fraction_used.py data/   (processes all files)
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


def compute_fraction_used(snapshots, N):
    """
    Compute F_u(t) = N_u(t) / N for each snapshot.
    """
    times = []
    fu = []
    
    for t, particles in snapshots:
        n_used = sum(1 for p in particles if p['state'] == 'U')
        times.append(t)
        fu.append(n_used / N)
    
    return np.array(times), np.array(fu)


def detect_steady_state(times, fu, threshold_fraction=0.05):
    """
    Detect when the system reaches steady state.
    
    Method: Find the first time when F_u stays within a band of width 
    `threshold_fraction * F_u_final` around the final average value.
    
    Returns:
        t_ss: time to reach steady state
        f_est: steady-state value of F_u
        f_est_std: standard deviation of F_u in steady state
    """
    if len(times) < 10 or np.max(fu) < 1e-10:
        return times[-1] if len(times) > 0 else 0.0, 0.0, 0.0
    
    # Steady-state value: mean of last 25%
    ss_start = int(0.75 * len(fu))
    f_est = np.mean(fu[ss_start:])
    f_est_std = np.std(fu[ss_start:])
    
    # Time to steady state: first time F_u enters the band and stays
    band = max(f_est_std * 2, f_est * threshold_fraction, 0.01)
    
    t_ss = times[-1]  # default to end of simulation
    window = max(1, len(fu) // 20)
    
    for i in range(len(fu) - window):
        segment = fu[i:i + window]
        if np.all(np.abs(segment - f_est) < band):
            t_ss = times[i]
            break
    
    return t_ss, f_est, f_est_std


def plot_fraction_used(metadata=None, snapshots=None, filepath=None,
                       output_dir="graphics/output"):
    """
    Plot F_u(t) for a single simulation file.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    if metadata is None or snapshots is None:
        if filepath is None:
            print("  [1.3] No data provided.")
            return
        metadata, snapshots, _ = parse_simulation_file(filepath)
    
    N = int(metadata['N'])
    times, fu = compute_fraction_used(snapshots, N)
    
    if len(times) == 0:
        print("  [1.3] No snapshots found.")
        return
    
    t_ss, f_est, f_est_std = detect_steady_state(times, fu)
    
    # Moving average for smoothing
    window = max(1, len(fu) // 20)
    if window > 1:
        fu_smooth = np.convolve(fu, np.ones(window) / window, mode='valid')
        t_smooth = times[:len(fu_smooth)]
    else:
        fu_smooth = fu
        t_smooth = times
    
    # ── Plot ──────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 7))
    
    ax.plot(times, fu, alpha=0.4, color='#7B1FA2', linewidth=0.8,
            label='$F_u(t)$ instantáneo')
    
    if window > 1:
        ax.plot(t_smooth, fu_smooth, color='#4A148C', linewidth=2.0,
                label=f'$F_u(t)$ promedio móvil (ventana={window})')
    
    # Steady-state band
    ax.axhline(y=f_est, color='#E91E63', linestyle='--', linewidth=1.5,
               label=f'$F_{{est}} = {f_est:.4f} \\pm {f_est_std:.4f}$')
    ax.axhspan(f_est - f_est_std, f_est + f_est_std, alpha=0.1, color='#E91E63')
    
    # Time to steady state
    ax.axvline(x=t_ss, color='#FF9800', linestyle=':', linewidth=1.5,
               label=f'$t_{{est}} = {t_ss:.3f}$ s')
    
    ax.set_xlabel('Tiempo $t$ [s]', fontsize=14)
    ax.set_ylabel('Fracción de partículas usadas $F_u(t) = N_u(t)/N$', fontsize=14)
    ax.set_title(f'Inciso 1.3: Evolución temporal de $F_u(t)$ ($N={N}$)',
                 fontsize=16, fontweight='bold')
    ax.legend(fontsize=11, loc='best')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', labelsize=12)
    ax.set_ylim(-0.05, max(1.05, np.max(fu) * 1.1))
    ax.set_xlim(times[0], times[-1])
    
    plt.tight_layout()
    outpath = os.path.join(output_dir, f"inciso_1_3_fraction_used_N{N}.png")
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [1.3] Saved: {outpath}")
    print(f"  [1.3] N={N}: F_est = {f_est:.4f} ± {f_est_std:.4f}, t_est = {t_ss:.3f} s")
    
    return t_ss, f_est, f_est_std


def plot_fest_vs_N(files_by_N, output_dir="graphics/output"):
    """
    Plot F_est(N) and t_est(N) — steady-state fraction and time as functions of N.
    Required by the enunciado: "Reportar tiempo al estacionario y valor del 
    estacionario alcanzado (Fest) en función de N."
    """
    os.makedirs(output_dir, exist_ok=True)
    
    N_values = sorted(files_by_N.keys())
    fest_means = []
    fest_stds = []
    test_means = []
    test_stds = []
    
    for N in N_values:
        fest_list = []
        test_list = []
        
        for entry in files_by_N[N]:
            fpath = entry[0]
            meta = entry[1]
            snaps = entry[2]
            
            n = int(meta['N'])
            times, fu = compute_fraction_used(snaps, n)
            t_ss, f_est, _ = detect_steady_state(times, fu)
            fest_list.append(f_est)
            test_list.append(t_ss)
        
        fest_means.append(np.mean(fest_list))
        fest_stds.append(np.std(fest_list))
        test_means.append(np.mean(test_list))
        test_stds.append(np.std(test_list))
    
    # ── Plot F_est(N) ─────────────────────────────────────────────────────
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))
    
    ax1.errorbar(N_values, fest_means, yerr=fest_stds, fmt='o-', capsize=5,
                 color='#7B1FA2', markerfacecolor='#CE93D8', markersize=8,
                 linewidth=2, label=r'$\langle F_{est} \rangle \pm \sigma$')
    ax1.set_xlabel('Número de partículas $N$', fontsize=14)
    ax1.set_ylabel(r'$F_{est}$ (fracción estacionaria)', fontsize=14)
    ax1.set_title(r'Inciso 1.3: $F_{est}$ vs $N$', fontsize=16, fontweight='bold')
    ax1.legend(fontsize=12)
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.tick_params(axis='both', labelsize=12)
    
    ax2.errorbar(N_values, test_means, yerr=test_stds, fmt='s-', capsize=5,
                 color='#FF6F00', markerfacecolor='#FFB74D', markersize=8,
                 linewidth=2, label=r'$\langle t_{est} \rangle \pm \sigma$')
    ax2.set_xlabel('Número de partículas $N$', fontsize=14)
    ax2.set_ylabel(r'Tiempo al estacionario $t_{est}$ [s]', fontsize=14)
    ax2.set_title(r'Inciso 1.3: $t_{est}$ vs $N$', fontsize=16, fontweight='bold')
    ax2.legend(fontsize=12)
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.tick_params(axis='both', labelsize=12)
    
    plt.tight_layout()
    outpath = os.path.join(output_dir, "inciso_1_3_fest_vs_N.png")
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [1.3] Saved: {outpath}")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python plot_fraction_used.py <sim_file_or_directory>")
        sys.exit(1)
    
    path = sys.argv[1]
    
    if os.path.isfile(path):
        metadata, snapshots, _ = parse_simulation_file(path)
        plot_fraction_used(metadata, snapshots, path)
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
        
        # Plot F_u(t) for each N (using first realization)
        for N in sorted(files_by_N.keys()):
            fpath, meta, snaps, _ = files_by_N[N][0]
            plot_fraction_used(meta, snaps, fpath)
        
        # Plot F_est(N) and t_est(N)
        if len(files_by_N) > 1:
            plot_fest_vs_N(files_by_N)
    else:
        print(f"Error: '{path}' not found.")
        sys.exit(1)
