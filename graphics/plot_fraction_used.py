"""
plot_fraction_used.py - Inciso 1.3
Generates: Temporal evolution of the fraction of used particles F_u(t).

Per the enunciado:
- "Estudiar la evolución temporal de la fracción de partículas usadas: Fu(t) = Nu(t)/N"
- "Reportar tiempo al estacionario y valor del estacionario alcanzado (Fest) en función de N"

This script generates:
1. F_u(t) vs t for each realization when a single file is provided
2. F_u(t) vs t for all realizations of a given N
3. A collapsed F_u(t) plot grouping realizations by N

Usage:
    python graphics/plot_fraction_used.py data/sim_300N_20260409_235938_s42.txt
    python graphics/plot_fraction_used.py data/   (processes all files)
"""

import sys
import os
import glob
import re
import numpy as np
import matplotlib.pyplot as plt

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


def describe_run(fpath, run_index):
    """Build a readable label for one realization."""
    basename = os.path.basename(fpath)
    seed_match = re.search(r'_s(\d+)\.txt$', basename)
    if seed_match:
        return f'Realizacion {run_index} (seed={seed_match.group(1)})'
    return f'Realizacion {run_index}'


def extract_fu_series(entry):
    """Return cached F_u(t) arrays for one realization."""
    return (
        np.array(entry["fu"]["times"], dtype=float),
        np.array(entry["fu"]["values"], dtype=float),
    )


def get_fraction_ylim(fu, fu_std=None):
    """
    Pick a y-axis range that matches the actual scale of F_u(t).
    """
    upper = float(np.max(fu))
    if fu_std is not None:
        upper = max(upper, float(np.max(fu + fu_std)))

    if upper < 1e-12:
        return 0.0, 0.05

    margin = max(upper * 0.15, 0.01)
    return 0.0, upper + margin


def plot_fraction_used(entry=None, filepath=None, output_dir="graphics/output"):
    """
    Plot F_u(t) for a single simulation file.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    if entry is None:
        if filepath is None:
            print("  [1.3] No data provided.")
            return
        entry = load_analysis_file(filepath)
    
    metadata = entry["metadata"]
    N = int(metadata['N'])
    times, fu = extract_fu_series(entry)
    
    if len(times) == 0:
        print("  [1.3] No snapshots found.")
        return
    
    # ── Plot ──────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 7))

    ax.step(times, fu, where='post', color='#5E35B1', linewidth=2.0,
            label='$F_u(t)$')
    
    ax.set_xlabel('Tiempo $t$ [s]', fontsize=14)
    ax.set_ylabel('Fracción de partículas usadas $F_u(t) = N_u(t)/N$', fontsize=14)
    ax.set_title(f'Inciso 1.3: Evolución temporal de $F_u(t)$ ($N={N}$)',
                 fontsize=16, fontweight='bold')
    ax.legend(fontsize=11, loc='upper left')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', labelsize=12)
    ymin, ymax = get_fraction_ylim(fu)
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(times[0], times[-1])
    
    plt.tight_layout()
    outpath = os.path.join(output_dir, f"inciso_1_3_fraction_used_N{N}.png")
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [1.3] Saved: {outpath}")
    print(f"  [1.3] N={N}: gráfico temporal de F_u(t) generado sin estimación automática del estacionario.")


def plot_fraction_used_realizations(entries, output_dir="graphics/output"):
    """
    Plot all F_u(t) realizations for a given N on the same axes.
    """
    os.makedirs(output_dir, exist_ok=True)

    metadata = entries[0]["metadata"]
    N = int(metadata['N'])
    n_runs = len(entries)
    colors = plt.cm.tab10(np.linspace(0, 1, max(n_runs, 1)))
    max_fu = 0.0
    t_max = 0.0

    fig, ax = plt.subplots(figsize=(13, 7))

    for run_index, entry in enumerate(entries, start=1):
        fpath = entry["source_path"]
        times, fu = extract_fu_series(entry)
        if len(times) == 0:
            continue
        max_fu = max(max_fu, float(np.max(fu)))
        t_max = max(t_max, float(times[-1]))
        ax.step(times, fu, where='post', alpha=0.85, linewidth=1.6,
                color=colors[run_index - 1],
                label=f'{describe_run(fpath, run_index)}: $F_u(t)$')

    ax.set_xlabel('Tiempo $t$ [s]', fontsize=14)
    ax.set_ylabel('Fracción de partículas usadas $F_u(t) = N_u(t)/N$', fontsize=14)
    ax.set_title(f'Inciso 1.3: $F_u(t)$ para $N={N}$\n'
                 f'{n_runs} realizaciones superpuestas',
                 fontsize=16, fontweight='bold')
    ax.legend(fontsize=10, loc='center left', bbox_to_anchor=(1.02, 0.5),
              borderaxespad=0.0)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', labelsize=12)
    ymin, ymax = get_fraction_ylim(np.array([0.0, max_fu]))
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(0.0, t_max)

    plt.tight_layout(rect=(0, 0, 0.88, 1))
    outpath = os.path.join(output_dir, f"inciso_1_3_fraction_used_N{N}.png")
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [1.3] Saved: {outpath}")
    print(f"  [1.3] N={N}: gráfico con {n_runs} realizaciones superpuestas.")


def plot_fraction_used_collapsed(files_by_N, output_dir="graphics/output"):
    """
    Plot all F_u(t) curves together, grouping realizations by N using color.
    """
    os.makedirs(output_dir, exist_ok=True)

    N_values = sorted(files_by_N.keys())
    colors = plt.cm.viridis(np.linspace(0, 1, max(len(N_values), 1)))
    max_fu = 0.0
    t_max = 0.0

    fig, ax = plt.subplots(figsize=(12, 8))

    for color, N in zip(colors, N_values):
        first_for_N = True
        for entry in files_by_N[N]:
            times, fu = extract_fu_series(entry)
            if len(times) == 0:
                continue
            max_fu = max(max_fu, float(np.max(fu)))
            t_max = max(t_max, float(times[-1]))
            label = f'$N={N}$' if first_for_N else None
            ax.step(times, fu, where='post', color=color, alpha=0.45,
                    linewidth=1.4, label=label)
            first_for_N = False

    ax.set_xlabel('Tiempo $t$ [s]', fontsize=14)
    ax.set_ylabel('Fracción de partículas usadas $F_u(t) = N_u(t)/N$', fontsize=14)
    ax.set_title('Inciso 1.3: $F_u(t)$ colapsado por $N$\n'
                 'cada color agrupa todas las realizaciones de un mismo N',
                 fontsize=16, fontweight='bold')
    ax.legend(fontsize=11, loc='upper left')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', labelsize=12)
    ymin, ymax = get_fraction_ylim(np.array([0.0, max_fu]))
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(0.0, t_max)

    plt.tight_layout()
    outpath = os.path.join(output_dir, "inciso_1_3_fraction_used_collapsed.png")
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [1.3] Saved: {outpath}")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python plot_fraction_used.py <sim_file_or_directory>")
        sys.exit(1)
    
    path = sys.argv[1]
    
    if os.path.isfile(path):
        plot_fraction_used(filepath=path)
    elif os.path.isdir(path):
        entries = load_analysis_entries(path)
        if not entries:
            print("No simulation files found.")
            sys.exit(1)
        files_by_N = group_entries_by_N(entries)
        
        # Plot F_u(t) realizations for each N and a collapsed summary across N.
        for N in sorted(files_by_N.keys()):
            plot_fraction_used_realizations(files_by_N[N])
        if len(files_by_N) > 1:
            plot_fraction_used_collapsed(files_by_N)
        
    else:
        print(f"Error: '{path}' not found.")
        sys.exit(1)
