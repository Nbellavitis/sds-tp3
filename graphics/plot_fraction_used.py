"""
plot_fraction_used.py - Inciso 1.3
Generates: Temporal evolution of the fraction of used particles F_u(t).

Per the enunciado:
- "Estudiar la evolución temporal de la fracción de partículas usadas: Fu(t) = Nu(t)/N"
- "Reportar tiempo al estacionario y valor del estacionario alcanzado (Fest) en función de N"

When a paired events_*.txt transition log exists, the cache builder reconstructs
F_u(t) from exact F->U / U->F transitions. Otherwise it falls back to the saved
snapshots.

This script generates:
1. F_u(t) vs t for one simulation when a single file is provided
2. F_u(t) vs t for all realizations of each N when a directory is provided
3. Summary plots for t_est(N) and F_est(N)

Usage:
    python graphics/plot_fraction_used.py data/sim_300N_20260409_235938_s42.txt
    python graphics/plot_fraction_used.py data/   (processes all files)
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from analysis_cache import group_entries_by_N, load_analysis_entries, load_analysis_file
from plot_style import apply_plot_style, get_distinct_series_styles


DISPLAY_PROFILE_N_VALUES = {100,150,200,250, 300,350,400,450, 500,550,600,650,700,750, 800}


apply_plot_style()


MANUAL_T_EST_BY_N = {
    50: 200.0,
    100: 200.0,
    150: 300.0,
    200: 300.0,
    250: 400.0,
    300: 450.0,
    350: 600.0,
    400: 650.0,
    450: 750.0,
    500: 850.0,
    550: 900.0,
    600: 1000.0,
    650: 1100.0,
    700: 1200.0,
    750: 1300.0,
    800: 1550.0,
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


def describe_run(_fpath, run_index):
    """Build a readable label for one realization."""
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


def normalize_entries(entries_or_entry):
    """Always return a list of analysis entries."""
    if entries_or_entry is None:
        return []
    if isinstance(entries_or_entry, dict):
        return [entries_or_entry]
    return list(entries_or_entry)


def compute_f_est_from_entries(entries_or_entry, t_est):
    """Average all sampled F_u values with t >= t_est."""
    if t_est is None:
        return None

    stationary_values = []
    for entry in normalize_entries(entries_or_entry):
        times, fu = extract_fu_series(entry)
        if len(times) == 0:
            continue
        stationary_values.extend(fu[times >= t_est].tolist())

    if not stationary_values:
        return None

    return float(np.mean(stationary_values))


def collect_stationary_values(entries_or_entry, t_est):
    """Return the pooled stationary-value bag for one N."""
    if t_est is None:
        return np.array([], dtype=float)

    stationary_values = []
    for entry in normalize_entries(entries_or_entry):
        times, fu = extract_fu_series(entry)
        if len(times) == 0:
            continue
        stationary_values.extend(fu[times >= t_est].tolist())
    return np.array(stationary_values, dtype=float)


def get_stationary_values(entries_or_entry, N):
    """Return manual t_est and data-driven F_est for one N."""
    t_est = MANUAL_T_EST_BY_N.get(N)
    f_est = compute_f_est_from_entries(entries_or_entry, t_est)
    return t_est, f_est


def add_stationary_reference_lines(ax, N, entries_or_entry):
    """Draw manual t_est and data-driven F_est references when available."""
    t_est, f_est = get_stationary_values(entries_or_entry, N)

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
    t_est, f_est = add_stationary_reference_lines(ax, N, entry)
    
    ax.set_xlabel('Tiempo $t$ [s]', fontsize=17)
    ax.set_ylabel('Fracción de partículas usadas $F_u(t) = N_u(t)/N$', fontsize=17)
    ax.legend(fontsize=13, loc='upper left')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', labelsize=14)
    ylim_samples = np.array([0.0, float(np.max(fu)), f_est if f_est is not None else 0.0])
    ymin, ymax, ystep = get_fraction_ylim(ylim_samples)
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(times[0], times[-1])
    ax.yaxis.set_major_locator(MultipleLocator(ystep))
    
    plt.tight_layout()
    outpath = os.path.join(output_dir, f"inciso_1_3_fraction_used_N{N}.png")
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [1.3] Saved: {outpath}")
    print(f"  [1.3] N={N}: gráfico temporal de F_u(t) con t_est={t_est} y F_est={f_est} calculado promediando F_u(t>=t_est).")


def plot_fraction_used_realizations(entries, output_dir="graphics/output"):
    """
    Plot all F_u(t) realizations for a given N on the same axes.
    """
    os.makedirs(output_dir, exist_ok=True)

    metadata = entries[0]["metadata"]
    N = int(metadata['N'])
    n_runs = len(entries)
    styles = get_distinct_series_styles(n_runs)
    max_fu = 0.0
    t_max = 0.0

    fig, ax = plt.subplots(figsize=(13, 7))

    for run_index, entry in enumerate(entries, start=1):
        times, fu = extract_fu_series(entry)
        if len(times) == 0:
            continue
        max_fu = max(max_fu, float(np.max(fu)))
        t_max = max(t_max, float(times[-1]))
        ax.step(times, fu, where='post', alpha=0.85, linewidth=1.6,
                color=styles[run_index - 1]["color"],
                label=describe_run(entry["source_path"], run_index))

    t_est, f_est = add_stationary_reference_lines(ax, N, entries)

    ax.set_xlabel('Tiempo $t$ [s]', fontsize=17)
    ax.set_ylabel('Fracción de partículas usadas $F_u(t)$', fontsize=17)
    ax.legend(fontsize=12, loc='upper right')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', labelsize=14)
    ylim_samples = np.array([0.0, max_fu, f_est if f_est is not None else 0.0])
    ymin, ymax, ystep = get_fraction_ylim(ylim_samples)
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(0.0, t_max)
    ax.yaxis.set_major_locator(MultipleLocator(ystep))

    plt.tight_layout()
    outpath = os.path.join(output_dir, f"inciso_1_3_fraction_used_N{N}.png")
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [1.3] Saved: {outpath}")
    print(f"  [1.3] N={N}: gráfico con {n_runs} realizaciones superpuestas y F_est={f_est} calculado como promedio de F_u(t>=t_est).")


def get_available_manual_N_values(files_by_N):
    """Return N values present in the data and with an available manual t_est."""
    return [
        N for N in sorted(files_by_N.keys())
        if N in MANUAL_T_EST_BY_N
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

    ax.set_xlabel('Número de partículas $N$', fontsize=17)
    ax.set_ylabel(r'Tiempo estacionario $t_{est}$', fontsize=17)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', labelsize=14)
    ax.set_xticks(N_values)
    ax.set_ylim(bottom=0.0)

    plt.tight_layout()
    outpath = os.path.join(output_dir, 'inciso_1_3_t_est_vs_N.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [1.3] Saved: {outpath}")


def plot_f_est_vs_N(files_by_N, output_dir="graphics/output"):
    """Plot data-driven F_est values as a function of N."""
    os.makedirs(output_dir, exist_ok=True)

    candidate_N_values = get_available_manual_N_values(files_by_N)
    if not candidate_N_values:
        print("  [1.3] No manual t_est values available for the loaded N set.")
        return

    N_values = []
    f_values = []
    f_stds = []
    for N in candidate_N_values:
        stationary_values = collect_stationary_values(files_by_N[N], MANUAL_T_EST_BY_N[N])
        if stationary_values.size == 0:
            continue
        f_est = float(np.mean(stationary_values))
        f_std = float(np.std(stationary_values))
        N_values.append(N)
        f_values.append(f_est)
        f_stds.append(f_std)

    if not N_values:
        print("  [1.3] No F_est values could be computed from the loaded data.")
        return

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.errorbar(
        N_values,
        f_values,
        yerr=f_stds,
        fmt='o-',
        capsize=5,
        linewidth=2.0,
        markersize=7,
        color='#D81B60',
        ecolor='#F8BBD0',
        markerfacecolor='#F48FB1',
        markeredgecolor='#880E4F',
    )

    ax.set_xlabel('Número de partículas $N$', fontsize=17)
    ax.set_ylabel(r'Fracción estacionaria media $\langle F_{est} \rangle$', fontsize=17)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', labelsize=14)
    ax.set_xticks(N_values)
    ax.set_ylim(bottom=0.0)

    plt.tight_layout()
    outpath = os.path.join(output_dir, 'inciso_1_3_fest_vs_N.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [1.3] Saved: {outpath}")
    print("  [1.3] Barras de error: desvío estándar poblacional de la bolsa estacionaria por N.")


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
        selected_profile_ns = [N for N in sorted(files_by_N.keys()) if N in DISPLAY_PROFILE_N_VALUES]

        for N in selected_profile_ns:
            plot_fraction_used_realizations(files_by_N[N])
        plot_t_est_vs_N(files_by_N)
        plot_f_est_vs_N(files_by_N)

    else:
        print(f"Error: '{path}' not found.")
        sys.exit(1)
