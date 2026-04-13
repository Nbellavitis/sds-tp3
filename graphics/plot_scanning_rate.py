"""
plot_scanning_rate.py - Inciso 1.2
Generates: Scanning rate <J>(N) with error bars.

The scanning rate J is obtained by linear interpolation of C_fc(t),
the cumulative count of fresh particles that contacted the central obstacle.
J = slope of the linear fit of C_fc(t).

C_fc(t) is read from explicit 'E' event lines in the output file (exact timestamps),
NOT from snapshot state diffing (which would miss fast F->U->F cycles).

The expensive extraction is delegated to a Java cache builder; Python only
loads the compact cached observables and plots them.

Usage:
    python graphics/plot_scanning_rate.py
    python graphics/plot_scanning_rate.py data/
"""

import sys
import os
import argparse
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analysis_cache import group_entries_by_N, load_analysis_entries, load_analysis_file


def parse_simulation_file(filepath):
    """
    Parse simulation file. Returns:
        metadata: dict
        snapshots: list of (time, particles)
        events: list of (time, particle_id) for F->U transitions
    """
    metadata = {}
    snapshots = []
    events = []  # F->U transition events with exact times
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    i = 0
    total = len(lines)
    
    # Parse header/metadata
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
            # Event line: E <time> <particle_id>
            parts = line.split()
            events.append((float(parts[1]), int(parts[2])))
        i += 1
    
    return metadata, snapshots, events


def compute_cfc_from_events(events, t_final):
    """
    Compute C_fc(t) from explicit F->U event lines.
    
    Returns:
        times: array of event times (including t=0 and t=t_final)
        cfc: cumulative count at each time
    """
    # Sort events by time
    events_sorted = sorted(events, key=lambda e: e[0])
    
    times = [0.0]
    cfc = [0]
    
    count = 0
    for t, pid in events_sorted:
        count += 1
        times.append(t)
        cfc.append(count)
    
    # Add final point
    if len(times) == 0 or times[-1] < t_final:
        times.append(t_final)
        cfc.append(count)
    
    return np.array(times), np.array(cfc)


def compute_cfc_from_snapshots(snapshots):
    """
    Fallback: Compute C_fc(t) from snapshot state diffing.
    Used when no 'E' event lines exist in the file.
    
    NOTE: This can MISS events where a particle goes F->U->F between snapshots.
    """
    times = []
    cfc_values = []
    cumulative = 0
    prev_states = None
    
    for t, particles in snapshots:
        current_states = {p['id']: p['state'] for p in particles}
        
        if prev_states is not None:
            for pid, state in current_states.items():
                if pid in prev_states and prev_states[pid] == 'F' and state == 'U':
                    cumulative += 1
        
        times.append(t)
        cfc_values.append(cumulative)
        prev_states = current_states
    
    return np.array(times), np.array(cfc_values)


def compute_scanning_rate(times, cfc):
    """
    Compute scanning rate J as the slope of the linear fit of C_fc(t).
    The enunciado says: "interpolar linealmente C_fc(t), su pendiente será J".
    We fit the ENTIRE curve (C_fc is cumulative so it should be roughly linear).
    
    Returns: J (slope), intercept, r_squared
    """
    if len(times) < 2:
        return 0.0, 0.0, 0.0
    
    # Linear fit over the full time range
    coeffs = np.polyfit(times, cfc, 1)
    J = coeffs[0]  # slope = scanning rate
    intercept = coeffs[1]
    
    # R-squared
    c_pred = np.polyval(coeffs, times)
    ss_res = np.sum((cfc - c_pred) ** 2)
    ss_tot = np.sum((cfc - np.mean(cfc)) ** 2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    
    return J, intercept, r_squared


def describe_run(fpath, run_index):
    """Build a human-readable label for a realization."""
    basename = os.path.basename(fpath)
    seed_match = re.search(r'_s(\d+)\.txt$', basename)
    if seed_match:
        return f'Realizacion {run_index} (seed={seed_match.group(1)})'
    return f'Realizacion {run_index}'


def plot_scanning_rate(files_by_N=None, data_dir="data", output_dir="graphics/output"):
    """
    Plot <J>(N) with error bars from multiple realizations.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    if files_by_N is None:
        entries = load_analysis_entries(data_dir)
        if not entries:
            print("  [1.2] No simulation files found.")
            return
        files_by_N = group_entries_by_N(entries)
    
    N_values = sorted(files_by_N.keys())
    J_means = []
    J_stds = []
    
    for N in N_values:
        J_list = []
        for entry in files_by_N[N]:
            J_list.append(float(entry["cfc"]["J"]))
        
        J_means.append(np.mean(J_list))
        J_stds.append(np.std(J_list))
        print(f"  [1.2] N={N}: <J> = {np.mean(J_list):.4f} ± {np.std(J_list):.4f} contacts/s "
              f"({len(J_list)} runs)")
    
    # ── Plot <J>(N) ──────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 7))
    
    ax.errorbar(N_values, J_means, yerr=J_stds, fmt='s-', capsize=5,
                color='#9C27B0', ecolor='#CE93D8', markerfacecolor='#6A1B9A',
                markeredgecolor='#4A148C', markersize=8, linewidth=2)
    
    ax.set_xlabel('Número de partículas $N$', fontsize=14)
    ax.set_ylabel(r'Scanning rate $\langle J \rangle$ [contactos/s]', fontsize=14)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', labelsize=12)
    ax.ticklabel_format(axis='y', style='scientific', scilimits=(0, 3))
    
    plt.tight_layout()
    outpath = os.path.join(output_dir, "inciso_1_2_scanning_rate.png")
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [1.2] Saved: {outpath}")
    
    # ── Also plot C_fc(t) for the largest N ──────────────────────────────
    max_N = max(N_values)
    fig2, ax2 = plt.subplots(figsize=(13, 7))
    
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(files_by_N[max_N]), 1)))

    for run_index, entry in enumerate(files_by_N[max_N], start=1):
        fpath = entry["source_path"]
        meta = entry["metadata"]
        t_final = float(meta.get('t_final', 5.0))
        times = np.array(entry["cfc"]["times"], dtype=float)
        cfc = np.array(entry["cfc"]["values"], dtype=float)
        J = float(entry["cfc"]["J"])
        intercept = float(entry["cfc"]["intercept"])
        color = colors[run_index - 1]

        ax2.step(times, cfc, where='post', alpha=0.8, linewidth=1.6, color=color)
        
        # Linear fit line
        t_fit = np.linspace(0, t_final, 100)
        ax2.plot(t_fit, J * t_fit + intercept, '--', alpha=0.8, color=color)
    
    ax2.set_xlabel('Tiempo $t$ [s]', fontsize=14)
    ax2.set_ylabel('$C_{fc}(t)$ [contactos acumulados]', fontsize=14)
    ax2.set_title(rf'$C_{{fc}}(t)$ para $N={max_N}$', fontsize=14)
    ax2.text(
        0.02,
        0.98,
        "Lineas solidas: contactos acumulados\nLineas punteadas: ajuste lineal",
        transform=ax2.transAxes,
        va='top',
        ha='left',
        fontsize=11,
        bbox=dict(boxstyle='round,pad=0.35', facecolor='white', alpha=0.9, edgecolor='#BDBDBD'),
    )
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.tick_params(axis='both', labelsize=12)
    
    plt.tight_layout()
    outpath2 = os.path.join(output_dir, "inciso_1_2_cfc_curve.png")
    plt.savefig(outpath2, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [1.2] Saved: {outpath2}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot scanning rate <J>(N)')
    parser.add_argument('path', nargs='?', default='data',
                        help='Directory with sim_*.txt files or a single simulation file')
    parser.add_argument('--output-dir', default='graphics/output',
                        help='Directory to save the plots')
    args = parser.parse_args()

    if os.path.isfile(args.path):
        entry = load_analysis_file(args.path)
        files_by_N = {int(entry["metadata"]["N"]): [entry]}
        plot_scanning_rate(files_by_N=files_by_N, output_dir=args.output_dir)
    else:
        plot_scanning_rate(data_dir=args.path, output_dir=args.output_dir)
