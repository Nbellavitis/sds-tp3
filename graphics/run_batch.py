"""
run_batch.py - Batch runner for multiple N values and realizations.
Runs the Java simulation for different N values and collects timing data.
Generates the timing.txt file needed for Inciso 1.1.

Usage:
    python graphics/run_batch.py
    python graphics/run_batch.py --n-values 50 100 150 200 250 300 --runs 5
    python graphics/run_batch.py --n-start 50 --n-stop 500 --n-step 50 --runs 5
    python graphics/run_batch.py --for-1-1 --runs 5
    python graphics/run_batch.py --for-1-1 --n-max 650 --runs 5
    python graphics/run_batch.py --n-values 50 100 150 200 250 300 --runs 5 --t-final 8.0
"""

import subprocess
import os
import re
import argparse


def run_batch(n_values, runs_per_n, seed_base=42, t_final=5.0):
    """
    Run the Java simulation for various N values and collect timing data.
    Writes timing.txt to data/ directory.
    """
    os.makedirs("data", exist_ok=True)
    timing_file = os.path.join("data", "timing.txt")

    with open(timing_file, 'w') as tf:
        tf.write("# N time_ms\n")
        tf.write(f"# runs_per_n={runs_per_n}, seed_base={seed_base}, t_final={t_final}\n")
    
    for N in n_values:
        for run in range(runs_per_n):
            seed = seed_base + run
            print(f"Running N={N}, run {run+1}/{runs_per_n} (seed={seed})...")
            
            cmd = [
                "java", "-cp", "engine/target/classes",
                "ar.edu.itba.sds.tp3.EventDrivenSimulation",
                "-N", str(N), "-seed", str(seed), "-t_final", str(t_final)
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=".")
            
            if result.returncode != 0:
                print(f"  ERROR: {result.stderr}")
                continue
            
            # Parse execution time
            time_ms = None
            for line in result.stdout.split('\n'):
                if 'Completed in' in line:
                    match = re.search(r'([\d.]+)\s*ms', line)
                    if match:
                        time_ms = float(match.group(1))
            
            if time_ms is not None:
                with open(timing_file, 'a') as tf:
                    tf.write(f"{N} {time_ms:.2f}\n")
                print(f"  Completed: {time_ms:.2f} ms")
            else:
                print(f"  WARNING: Could not parse timing from output")
                print(f"  stdout: {result.stdout[:200]}")
    
    print(f"\nTiming data saved to: {timing_file}")


def build_n_values(args):
    """Resolve the list of N values from explicit values or a numeric range."""
    if args.for_1_1:
        if args.n_max < 50:
            raise ValueError("--n-max must be at least 50 for --for-1-1")
        return list(range(50, args.n_max + 1, 50))

    if args.n_values is not None:
        return args.n_values

    if args.n_step <= 0:
        raise ValueError("--n-step must be positive")
    if args.n_stop < args.n_start:
        raise ValueError("--n-stop must be greater than or equal to --n-start")

    return list(range(args.n_start, args.n_stop + 1, args.n_step))


def main():
    parser = argparse.ArgumentParser(description='Batch runner for simulations')
    parser.add_argument('--for-1-1', action='store_true',
                        help='Use the inciso 1.1 setup: N=50..Nmax with step 50 and t_f=5 s')
    parser.add_argument('--n-max', type=int, default=500,
                        help='Maximum N when using --for-1-1 (default: 500)')
    parser.add_argument('--n-values', nargs='+', type=int,
                        default=None,
                        help='Explicit N values to simulate')
    parser.add_argument('--n-start', type=int, default=50,
                        help='Start of the N range when --n-values is not provided')
    parser.add_argument('--n-stop', type=int, default=300,
                        help='End of the N range when --n-values is not provided')
    parser.add_argument('--n-step', type=int, default=50,
                        help='Step of the N range when --n-values is not provided')
    parser.add_argument('--runs', type=int, default=5,
                        help='Number of runs per N value')
    parser.add_argument('--seed', type=int, default=42,
                        help='Base seed for reproducibility')
    parser.add_argument('--t-final', type=float, default=5.0,
                        help='Simulation final time t_f in seconds')
    args = parser.parse_args()

    t_final = 5.0 if args.for_1_1 else args.t_final
    run_batch(build_n_values(args), args.runs, args.seed, t_final)


if __name__ == '__main__':
    main()
