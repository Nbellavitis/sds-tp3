"""
plot_fraction_used.py - Inciso 1.3
Generates: Temporal evolution of the fraction of used particles F_u(t).

Per the enunciado:
- "Estudiar la evolución temporal de la fracción de partículas usadas: Fu(t) = Nu(t)/N"
- "Reportar tiempo al estacionario y valor del estacionario alcanzado (Fest) en función de N"

This script generates:
1. F_u(t) vs t for one simulation when a single file is provided
2. A single overlay figure with all realizations across all N when a directory is provided

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


DEFAULT_T_STATIONARY = 400.0


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


def get_zoom_ylim(series_max):
    """Y-range tuned for visual plateau detection in F_u(t)."""
    upper = max(float(series_max), 1e-6)
    target = upper * 1.12
    if target <= 0.1:
        step = 0.01
    elif target <= 0.2:
        step = 0.02
    elif target <= 0.5:
        step = 0.05
    else:
        step = 0.1
    ymax = min(1.0, np.ceil(target / step) * step)
    return 0.0, float(ymax), float(step)


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
    ymin, ymax, ystep = get_fraction_ylim(fu)
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(times[0], times[-1])
    ax.yaxis.set_major_locator(MultipleLocator(ystep))
    
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
    ymin, ymax, ystep = get_fraction_ylim(np.array([0.0, max_fu]))
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(0.0, t_max)
    ax.yaxis.set_major_locator(MultipleLocator(ystep))

    plt.tight_layout(rect=(0, 0, 0.88, 1))
    outpath = os.path.join(output_dir, f"inciso_1_3_fraction_used_N{N}.png")
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [1.3] Saved: {outpath}")
    print(f"  [1.3] N={N}: gráfico con {n_runs} realizaciones superpuestas.")


def plot_fraction_used_across_N(files_by_N, output_dir="graphics/output"):
    """
    Plot F_u(t) superposing all N using raw realizations only.

    This view is intended to estimate the transient-to-stationary time by eye
    without averaging realizations beforehand.
    """
    os.makedirs(output_dir, exist_ok=True)

    N_values = sorted(files_by_N.keys())
    if not N_values:
        print("  [1.3] No N groups available for cross-N overlay.")
        return

    all_tmax = []
    all_series_max = []

    grouped_series = {}
    for N in N_values:
        runs = []
        for entry in files_by_N[N]:
            times, fu = extract_fu_series(entry)
            if len(times) == 0:
                continue
            runs.append((times, fu))
            all_tmax.append(float(times[-1]))
            all_series_max.append(float(np.max(fu)))
        grouped_series[N] = runs

    if not all_tmax:
        print("  [1.3] Could not build cross-N overlay (no time series found).")
        return

    t_max = max(all_tmax)

    fig, ax = plt.subplots(figsize=(12, 7))
    colors = plt.cm.viridis(np.linspace(0.08, 0.92, len(N_values)))

    legend_labels = []
    for idx, N in enumerate(N_values):
        runs = grouped_series[N]
        if not runs:
            continue

        color = colors[idx]
        for run_index, (times, fu) in enumerate(runs, start=1):
            label = None
            if run_index == 1:
                label = f'N={N} ({len(runs)} realizaciones)'
                legend_labels.append(label)
            ax.step(times, fu, where='post', linewidth=1.4, alpha=0.35,
                    color=color, label=label)

    ymin, ymax, ystep = get_zoom_ylim(max(all_series_max))
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(0.0, t_max)
    ax.yaxis.set_major_locator(MultipleLocator(ystep))

    ax.set_xlabel('Tiempo $t$ [s]', fontsize=14)
    ax.set_ylabel('Fracción de partículas usadas $F_u(t)$', fontsize=14)
    ax.set_title('Inciso 1.3: $F_u(t)$ superpuesto para distintos $N$\n'
                 'series crudas por realización (sin promedio previo, escala Y ampliada)',
                 fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', labelsize=12)
    ax.legend(fontsize=10, loc='center left', bbox_to_anchor=(1.02, 0.5),
              borderaxespad=0.0)

    plt.tight_layout(rect=(0, 0, 0.84, 1))
    outpath = os.path.join(output_dir, 'inciso_1_3_fraction_used_all_N_overlay.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [1.3] Saved: {outpath}")
    print(f"  [1.3] Overlay por N listo (sin promedio previo; Y in [0, {ymax:.3f}] para inspección visual del estacionario).")


def seed_from_entry(entry):
    """Extract seed from source filename when available."""
    src = os.path.basename(entry.get("source_path", ""))
    match = re.search(r'_s(\d+)\.txt$', src)
    return int(match.group(1)) if match else None


def pick_one_entry_per_N(files_by_N):
    """Pick one deterministic realization per N (lowest seed if present)."""
    selected = {}
    for N, entries in files_by_N.items():
        if not entries:
            continue
        selected_entry = min(
            entries,
            key=lambda e: (
                seed_from_entry(e) is None,
                seed_from_entry(e) if seed_from_entry(e) is not None else 10**18,
                e.get("source_path", ""),
            ),
        )
        selected[N] = selected_entry
    return selected


def plot_fraction_used_one_run_per_N(files_by_N, output_dir="graphics/output"):
    """
    Plot one raw realization per N on the same figure (no averaging).
    """
    os.makedirs(output_dir, exist_ok=True)

    selected = pick_one_entry_per_N(files_by_N)
    N_values = sorted(selected.keys())
    if not N_values:
        print("  [1.3] No entries available for one-run-per-N plot.")
        return

    colors = plt.cm.viridis(np.linspace(0.08, 0.92, len(N_values)))
    fig, ax = plt.subplots(figsize=(12, 7))

    max_fu = 0.0
    t_max = 0.0
    for idx, N in enumerate(N_values):
        entry = selected[N]
        times, fu = extract_fu_series(entry)
        if len(times) == 0:
            continue

        max_fu = max(max_fu, float(np.max(fu)))
        t_max = max(t_max, float(times[-1]))
        seed = seed_from_entry(entry)
        seed_label = f", seed={seed}" if seed is not None else ""
        ax.step(times, fu, where='post', linewidth=1.8, alpha=0.95,
                color=colors[idx], label=f'N={N}{seed_label}')

    ymin, ymax, ystep = get_zoom_ylim(max_fu)
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(0.0, t_max)
    ax.yaxis.set_major_locator(MultipleLocator(ystep))
    ax.set_xlabel('Tiempo $t$ [s]', fontsize=14)
    ax.set_ylabel('Fracción de partículas usadas $F_u(t)$', fontsize=14)
    ax.set_title('Inciso 1.3: una realización cruda por cada $N$\n'
                 'superposición para estimar visualmente el estacionario',
                 fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', labelsize=12)
    ax.legend(fontsize=10, loc='center left', bbox_to_anchor=(1.02, 0.5),
              borderaxespad=0.0)

    plt.tight_layout(rect=(0, 0, 0.84, 1))
    outpath = os.path.join(output_dir, 'inciso_1_3_fraction_used_one_run_per_N.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [1.3] Saved: {outpath}")
    print("  [1.3] Gráfico crudo por N listo (una realización por N, sin promedio).")


def compute_fest_for_run(times, fu, t_stationary=DEFAULT_T_STATIONARY):
    """Compute time-weighted stationary mean of F_u(t) over [t_stationary, t_final]."""
    if len(times) < 2:
        return np.nan

    t_start = max(float(t_stationary), float(times[0]))
    t_end = float(times[-1])
    if t_end <= t_start:
        return np.nan

    weighted_sum = 0.0
    total_dt = 0.0
    for i in range(len(times) - 1):
        seg_start = float(times[i])
        seg_end = float(times[i + 1])
        overlap_start = max(seg_start, t_start)
        overlap_end = min(seg_end, t_end)
        dt = overlap_end - overlap_start
        if dt > 0:
            weighted_sum += dt * float(fu[i])
            total_dt += dt

    if total_dt <= 0:
        return np.nan
    return weighted_sum / total_dt


def plot_fest_vs_N(files_by_N, output_dir="graphics/output", t_stationary=DEFAULT_T_STATIONARY):
    """
    Plot final result for inciso 1.3: <F_est> vs N with ensemble std error bars.
    """
    os.makedirs(output_dir, exist_ok=True)

    N_values = []
    fest_means = []
    fest_stds = []

    for N in sorted(files_by_N.keys()):
        fest_runs = []
        for entry in files_by_N[N]:
            times, fu = extract_fu_series(entry)
            fest_run = compute_fest_for_run(times, fu, t_stationary)
            if np.isfinite(fest_run):
                fest_runs.append(fest_run)

        if not fest_runs:
            print(f"  [1.3] WARNING: N={N} has no valid samples for t >= {t_stationary:g} s")
            continue

        fest_runs = np.array(fest_runs, dtype=float)
        N_values.append(N)
        fest_means.append(float(np.mean(fest_runs)))
        if len(fest_runs) > 1:
            fest_stds.append(float(np.std(fest_runs, ddof=1)))
        else:
            fest_stds.append(0.0)

        print(f"  [1.3] N={N}: <F_est> = {fest_means[-1]:.4f} ± {fest_stds[-1]:.4f} ({len(fest_runs)} runs)")

    if not N_values:
        print("  [1.3] No valid data to plot <F_est> vs N.")
        return

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.errorbar(
        N_values,
        fest_means,
        yerr=fest_stds,
        fmt='o-',
        linewidth=1.4,
        markersize=7,
        capsize=5,
        color='#6A1B9A',
        ecolor='#CE93D8',
        markerfacecolor='#8E24AA',
        markeredgecolor='#4A148C',
        label=r'$\langle F_{est} \rangle \pm \sigma$'
    )

    ax.set_xlabel('Número de partículas $N$', fontsize=14)
    ax.set_ylabel(r'Fracción estacionaria media $\langle F_{est} \rangle$', fontsize=14)
    ax.set_title(rf'Inciso 1.3: Resultado final $\langle F_{{est}} \rangle$ vs $N$ '
                 rf'($t_{{est}}={t_stationary:g}$ s)', fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', labelsize=12)
    ax.set_ylim(bottom=0.0)
    ax.legend(fontsize=11)

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

        # Directory mode: one raw realization per N in a single figure.
        plot_fraction_used_one_run_per_N(files_by_N)
        plot_fest_vs_N(files_by_N, t_stationary=DEFAULT_T_STATIONARY)

    else:
        print(f"Error: '{path}' not found.")
        sys.exit(1)
