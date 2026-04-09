"""
plot_execution_time.py - Inciso 1.1
Generates: Execution Time vs N (up to t_f = 5 s)

Reads multiple simulation output files from data/ directory,
extracts the wall-clock execution time from each, and plots
execution time vs number of particles N.

Usage:
    python graphics/plot_execution_time.py data/
"""

import sys
import os
import re
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def run_simulations_and_measure(N_values, runs_per_N=5, seed_base=42):
    """
    Run simulations for various N values and measure execution time.
    Returns dict: N -> list of execution times [ms]
    
    This function is called when generating the plot from scratch by running
    the Java simulation multiple times.
    """
    times_by_N = {}
    
    for N in N_values:
        times_by_N[N] = []
        for run in range(runs_per_N):
            seed = seed_base + run
            cmd = [
                "java", "-cp", "engine/target/classes",
                "ar.edu.itba.sds.tp3.EventDrivenSimulation",
                "-N", str(N), "-seed", str(seed)
            ]
            print(f"  Running N={N}, run {run+1}/{runs_per_N}...")
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=".")
            
            # Parse execution time from output
            for line in result.stdout.split('\n'):
                if 'Completed in' in line:
                    match = re.search(r'(\d+\.?\d*)\s*ms', line)
                    if match:
                        times_by_N[N].append(float(match.group(1)))
    
    return times_by_N


def parse_time_from_files(data_dir):
    """
    Alternative: parse execution times from a pre-generated timing file.
    Expected format: timing.txt with lines "N time_ms"
    """
    timing_file = os.path.join(data_dir, "timing.txt")
    times_by_N = {}
    
    if os.path.exists(timing_file):
        with open(timing_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    parts = line.split()
                    N = int(parts[0])
                    t = float(parts[1])
                    if N not in times_by_N:
                        times_by_N[N] = []
                    times_by_N[N].append(t)
    
    return times_by_N


def plot_execution_time(files_by_N=None, output_dir="graphics/output"):
    """
    Plot execution time vs N.
    
    Args:
        files_by_N: dict {N: [(filepath, metadata, snapshots), ...]}
                    If None, tries to read from timing.txt
        output_dir: directory to save the plot
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Try to read timing data
    times_by_N = parse_time_from_files("data")
    
    if not times_by_N:
        print("  [1.1] No timing data found. Run simulations first or provide timing.txt")
        print("  [1.1] Format of timing.txt: one line per run with 'N time_ms'")
        return
    
    N_values = sorted(times_by_N.keys())
    means = [np.mean(times_by_N[n]) for n in N_values]
    stds = [np.std(times_by_N[n]) for n in N_values]
    
    # ── Plot ──────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 7))
    
    ax.errorbar(N_values, means, yerr=stds, fmt='o-', capsize=5,
                color='#2196F3', ecolor='#90CAF9', markerfacecolor='#1565C0',
                markeredgecolor='#0D47A1', markersize=8, linewidth=2,
                label='Tiempo medio ± σ')
    
    ax.set_xlabel('Número de partículas $N$', fontsize=14)
    ax.set_ylabel('Tiempo de ejecución [ms]', fontsize=14)
    ax.set_title('Inciso 1.1: Tiempo de ejecución vs $N$ ($t_f = 5$ s)',
                 fontsize=16, fontweight='bold')
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', labelsize=12)
    
    # Use scientific notation on y-axis if values are large
    ax.ticklabel_format(axis='y', style='scientific', scilimits=(0, 3))
    
    plt.tight_layout()
    outpath = os.path.join(output_dir, "inciso_1_1_execution_time.png")
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [1.1] Saved: {outpath}")


if __name__ == '__main__':
    if len(sys.argv) > 1:
        plot_execution_time(output_dir=sys.argv[1] if len(sys.argv) > 1 else "graphics/output")
    else:
        plot_execution_time()
