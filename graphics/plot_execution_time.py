"""
plot_execution_time.py - Inciso 1.1
Generates: Execution Time vs N.

Reads multiple simulation output files from data/ directory,
extracts the wall-clock execution time from each, and plots
execution time vs number of particles N.

Usage:
    python graphics/plot_execution_time.py data/
"""

import os
import re
import subprocess
import argparse
import numpy as np
import matplotlib.pyplot as plt


REQUIRED_T_FINAL_1_1 = 5.0


def fit_exponential_model(x, y):
    """
    Fit y = A * exp(b x) using linear regression on log(y).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    positive = y > 0
    if np.count_nonzero(positive) < 2:
        return None

    log_y = np.log(y[positive])
    coeffs = np.polyfit(x[positive], log_y, 1)
    b = coeffs[0]
    log_a = coeffs[1]
    y_fit = np.exp(log_a + b * x)
    ss_res = np.sum((log_y - (log_a + b * x[positive])) ** 2)
    ss_tot = np.sum((log_y - np.mean(log_y)) ** 2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    return {
        "A": float(np.exp(log_a)),
        "b": float(b),
        "y_fit": y_fit,
        "r_squared_log": float(r_squared),
    }


def fit_power_law_model(x, y):
    """
    Fit y = C * x^alpha using linear regression on log(x), log(y).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    positive = (x > 0) & (y > 0)
    if np.count_nonzero(positive) < 2:
        return None

    log_x = np.log(x[positive])
    log_y = np.log(y[positive])
    coeffs = np.polyfit(log_x, log_y, 1)
    alpha = coeffs[0]
    log_c = coeffs[1]
    y_fit = np.exp(log_c) * np.power(x, alpha)
    ss_res = np.sum((log_y - (log_c + alpha * log_x)) ** 2)
    ss_tot = np.sum((log_y - np.mean(log_y)) ** 2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    return {
        "C": float(np.exp(log_c)),
        "alpha": float(alpha),
        "y_fit": y_fit,
        "r_squared_loglog": float(r_squared),
    }


def run_simulations_and_measure(N_values, runs_per_N=5, seed_base=42, t_final=5.0):
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
                "-N", str(N), "-seed", str(seed), "-t_final", str(t_final)
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


def parse_time_from_file(timing_file):
    """Parse one timing file with optional metadata header."""
    times_by_N = {}
    metadata = {}

    if not os.path.exists(timing_file):
        return times_by_N, metadata

    with open(timing_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                for key, value in re.findall(r'([A-Za-z_]+)=([^\s,]+)', line):
                    try:
                        metadata[key] = float(value)
                    except ValueError:
                        metadata[key] = value
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            N = int(parts[0])
            t = float(parts[1])
            times_by_N.setdefault(N, []).append(t)

    return times_by_N, metadata


def parse_time_from_files(data_dir):
    """
    Alternative: parse execution times from a pre-generated timing file.
    Expected format: timing.txt with lines "N time_ms"
    """
    preferred = os.path.join(data_dir, "timing_1_1.txt")
    fallback = os.path.join(data_dir, "timing.txt")

    if os.path.exists(preferred):
        times_by_N, metadata = parse_time_from_file(preferred)
        metadata["_timing_file"] = preferred
        return times_by_N, metadata

    times_by_N, metadata = parse_time_from_file(fallback)
    metadata["_timing_file"] = fallback
    return times_by_N, metadata


def plot_execution_time(files_by_N=None, data_dir="data", output_dir="graphics/output"):
    """
    Plot execution time vs N.
    
    Args:
        files_by_N: dict {N: [(filepath, metadata, snapshots), ...]}
                    If None, tries to read from timing.txt
        output_dir: directory to save the plot
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Try to read timing data
    times_by_N, timing_metadata = parse_time_from_files(data_dir)
    
    timing_file = timing_metadata.get('_timing_file', os.path.join(data_dir, 'timing.txt'))

    if not times_by_N:
        print("  [1.1] No timing data found. Run simulations first or provide timing.txt")
        print("  [1.1] Format of timing.txt: one line per run with 'N time_ms'")
        return

    t_final = timing_metadata.get('t_final')
    if not isinstance(t_final, (int, float)):
        print(f"  [1.1] WARNING: {timing_file} has no t_final metadata. "
              f"Inciso 1.1 requires t_f = {REQUIRED_T_FINAL_1_1:g} s.")
        print("  [1.1] Re-generate timing with: python graphics/run_batch.py --for-1-1 --runs 5")
        return

    if abs(float(t_final) - REQUIRED_T_FINAL_1_1) > 1e-9:
        print(f"  [1.1] WARNING: {timing_file} was generated with t_f={t_final:g} s. "
              f"Inciso 1.1 requires t_f = {REQUIRED_T_FINAL_1_1:g} s.")
        print("  [1.1] Re-generate timing with: python graphics/run_batch.py --for-1-1 --runs 5")
        return

    N_values = sorted(times_by_N.keys())
    means = np.array([np.mean(times_by_N[n]) for n in N_values], dtype=float)
    stds = np.array([np.std(times_by_N[n]) for n in N_values], dtype=float)
    N_array = np.array(N_values, dtype=float)
    exp_fit = fit_exponential_model(N_array, means)
    power_fit = fit_power_law_model(N_array, means)
    
    # ── Plot ──────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 7))
    
    ax.errorbar(N_values, means, yerr=stds, fmt='o-', capsize=5,
                color='#2196F3', ecolor='#90CAF9', markerfacecolor='#1565C0',
                markeredgecolor='#0D47A1', markersize=8, linewidth=2,
                label='Tiempo medio ± σ')
    if exp_fit is not None:
        ax.plot(N_array, exp_fit["y_fit"], '--', color='#D81B60', linewidth=2,
                label=(r'Ajuste exp.: $T(N)\approx %.2f\,e^{%.4fN}$'
                       % (exp_fit["A"], exp_fit["b"])))

    ax.set_xlabel('Número de partículas $N$', fontsize=14)
    ax.set_ylabel('Tiempo de ejecución [ms] (escala log)', fontsize=14)
    ax.set_yscale('log')
    title = f'Inciso 1.1: Tiempo de ejecución vs $N$ ($t_f = {float(t_final):g}$ s)'
    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', labelsize=12)
    plt.tight_layout()
    outpath = os.path.join(output_dir, "inciso_1_1_execution_time.png")
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [1.1] Saved: {outpath}")

    if exp_fit is not None:
        print(f"  [1.1] Exponential fit: T(N) ≈ {exp_fit['A']:.3f} * exp({exp_fit['b']:.5f} N)"
              f" with R²(log-space) = {exp_fit['r_squared_log']:.4f}")
    if power_fit is not None:
        print(f"  [1.1] Power-law fit: T(N) ≈ {power_fit['C']:.5f} * N^{power_fit['alpha']:.3f}"
              f" with R²(log-log) = {power_fit['r_squared_loglog']:.4f}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot execution time vs N')
    parser.add_argument('data_dir', nargs='?', default='data',
                        help='Directory containing timing.txt (default: data)')
    parser.add_argument('--output-dir', default='graphics/output',
                        help='Directory to save the plot')
    args = parser.parse_args()

    plot_execution_time(data_dir=args.data_dir, output_dir=args.output_dir)
