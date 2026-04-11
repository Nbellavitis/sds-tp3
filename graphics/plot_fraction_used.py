"""
plot_fraction_used.py - Inciso 1.3
Generates: Temporal evolution of the fraction of used particles F_u(t).

Per the enunciado:
- "Estudiar la evolución temporal de la fracción de partículas usadas: Fu(t) = Nu(t)/N"
- "Reportar tiempo al estacionario y valor del estacionario alcanzado (Fest) en función de N"

This script generates:
1. F_u(t) vs t for one simulation when a single file is provided
2. F_u(t) vs t for all realizations of each N when a directory is provided
3. Manual summary plots for t_est(N) and F_est(N)

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
from matplotlib.ticker import MultipleLocator

from analysis_cache import group_entries_by_N, load_analysis_entries, load_analysis_file


MANUAL_T_EST_BY_N = {
    50: 200.0,
    100: 200.0,
    150: 300.0,
    200: 300.0,
    250: 400.0,
    300: 400.0,
}

MANUAL_F_EST_BY_N = {
    50: 0.06,
    100: 0.08,
    150: 0.09,
    200: 0.11,
    250: 0.12,
    300: 0.13,
}


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
    Pick a y-axis range and tick spacing that avoid over-stretching
    very small F_u(t) values.
    """
    upper = float(np.max(fu))
    if fu_std is not None:
        upper = max(upper, float(np.max(fu + fu_std)))

    if upper < 1e-12:
        return 0.0, 0.05, 0.01

    target_step = max(upper, 0.04) / 4.0
    exponent = np.floor(np.log10(target_step))
    normalized = target_step / (10 ** exponent)
    if normalized <= 1.0:
        step = 1.0
    elif normalized <= 2.0:
        step = 2.0
    elif normalized <= 5.0:
        step = 5.0
    else:
        step = 10.0
    step *= 10 ** exponent

    ymax = max(5.0 * step, np.ceil(upper / step) * step)
    return 0.0, float(ymax), float(step)


def get_manual_stationary_values(N):
    """Return manually assigned stationary values for one N."""
    return MANUAL_T_EST_BY_N.get(N), MANUAL_F_EST_BY_N.get(N)


def add_stationary_reference_lines(ax, N):
    """Draw manual t_est and F_est references when available."""
    t_est, f_est = get_manual_stationary_values(N)

    if f_est is not None:
        ax.axhline(
            f_est,
            color="#D81B60",
            linestyle="--",
            linewidth=1.8,
            label=rf"$F_{{est}}={f_est:.2f}$",
        )

    if t_est is not None:
        ax.axvline(
            t_est,
            color="#FB8C00",
            linestyle="--",
            linewidth=1.8,
            label=rf"$t_{{est}}={t_est:g}$",
        )

    return t_est, f_est


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
    t_est, f_est = add_stationary_reference_lines(ax, N)
    
    ax.set_xlabel('Tiempo $t$ [s]', fontsize=14)
    ax.set_ylabel('Fracción de partículas usadas $F_u(t) = N_u(t)/N$', fontsize=14)
    ax.set_title(f'Inciso 1.3: Evolución temporal de $F_u(t)$ ($N={N}$)',
                 fontsize=16, fontweight='bold')
    ax.legend(fontsize=11, loc='upper left')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', labelsize=12)
    ylim_samples = np.array([0.0, float(np.max(fu)), f_est or 0.0])
    ymin, ymax, ystep = get_fraction_ylim(ylim_samples)
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(times[0], times[-1])
    ax.yaxis.set_major_locator(MultipleLocator(ystep))
    
    plt.tight_layout()
    outpath = os.path.join(output_dir, f"inciso_1_3_fraction_used_N{N}.png")
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [1.3] Saved: {outpath}")
    print(f"  [1.3] N={N}: gráfico temporal de F_u(t) con referencias manuales t_est={t_est} y F_est={f_est}.")


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

    t_est, f_est = add_stationary_reference_lines(ax, N)

    ax.set_xlabel('Tiempo $t$ [s]', fontsize=14)
    ax.set_ylabel('Fracción de partículas usadas $F_u(t) = N_u(t)/N$', fontsize=14)
    ax.set_title(f'Inciso 1.3: $F_u(t)$ para $N={N}$\n'
                 f'{n_runs} realizaciones superpuestas',
                 fontsize=16, fontweight='bold')
    ax.legend(fontsize=10, loc='center left', bbox_to_anchor=(1.02, 0.5),
              borderaxespad=0.0)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', labelsize=12)
    ylim_samples = np.array([0.0, max_fu, f_est or 0.0])
    ymin, ymax, ystep = get_fraction_ylim(ylim_samples)
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(0.0, t_max)
    ax.yaxis.set_major_locator(MultipleLocator(ystep))

    plt.tight_layout(rect=(0, 0, 0.88, 1))
    outpath = os.path.join(output_dir, f"inciso_1_3_fraction_used_N{N}.png")
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [1.3] Saved: {outpath}")
    print(f"  [1.3] N={N}: gráfico con {n_runs} realizaciones superpuestas y referencias t_est={t_est}, F_est={f_est}.")


def get_available_manual_N_values(files_by_N):
    """Return N values present in the data and in both manual estimate tables."""
    return [
        N for N in sorted(files_by_N.keys())
        if N in MANUAL_T_EST_BY_N and N in MANUAL_F_EST_BY_N
    ]


def plot_t_est_vs_N(files_by_N, output_dir="graphics/output"):
    """Plot manual t_est values as a function of N."""
    os.makedirs(output_dir, exist_ok=True)

    N_values = get_available_manual_N_values(files_by_N)
    if not N_values:
        print("  [1.3] No manual t_est values available for the loaded N set.")
        return

    t_values = [MANUAL_T_EST_BY_N[N] for N in N_values]

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(
        N_values,
        t_values,
        'o-',
        linewidth=2.0,
        markersize=7,
        color='#FB8C00',
        markerfacecolor='#FFB74D',
        markeredgecolor='#E65100',
    )

    ax.set_xlabel('Número de partículas $N$', fontsize=14)
    ax.set_ylabel(r'Tiempo estacionario $t_{est}$', fontsize=14)
    ax.set_title(r'Inciso 1.3: $t_{est}$ en función de $N$',
                 fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', labelsize=12)
    ax.set_xticks(N_values)
    ax.set_ylim(bottom=0.0)

    plt.tight_layout()
    outpath = os.path.join(output_dir, 'inciso_1_3_t_est_vs_N.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [1.3] Saved: {outpath}")


def plot_f_est_vs_N(files_by_N, output_dir="graphics/output"):
    """Plot manual F_est values as a function of N."""
    os.makedirs(output_dir, exist_ok=True)

    N_values = get_available_manual_N_values(files_by_N)
    if not N_values:
        print("  [1.3] No manual F_est values available for the loaded N set.")
        return

    f_values = [MANUAL_F_EST_BY_N[N] for N in N_values]

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(
        N_values,
        f_values,
        'o-',
        linewidth=2.0,
        markersize=7,
        color='#D81B60',
        markerfacecolor='#F48FB1',
        markeredgecolor='#880E4F',
    )

    ax.set_xlabel('Número de partículas $N$', fontsize=14)
    ax.set_ylabel(r'Valor estacionario $F_{est}$', fontsize=14)
    ax.set_title(r'Inciso 1.3: $F_{est}$ en función de $N$',
                 fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', labelsize=12)
    ax.set_xticks(N_values)
    ax.set_ylim(bottom=0.0)

    plt.tight_layout()
    outpath = os.path.join(output_dir, 'inciso_1_3_fest_vs_N.png')
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

        for N in sorted(files_by_N.keys()):
            plot_fraction_used_realizations(files_by_N[N])
        plot_t_est_vs_N(files_by_N)
        plot_f_est_vs_N(files_by_N)

    else:
        print(f"Error: '{path}' not found.")
        sys.exit(1)
