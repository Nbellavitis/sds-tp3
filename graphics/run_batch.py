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
    python graphics/run_batch.py --for-1-1 --n-max 650 --runs 5 --t-final 8.0
    python graphics/run_batch.py --n-values 50 100 150 200 250 300 --runs 5 --t-final 8.0
"""

import subprocess
import os
import re
import argparse


DEFAULT_TIMING_WINDOW_1_1 = 5.0


def parse_runtime_report(stdout):
    """Extract full-run timing and the dedicated inciso 1.1 timing from Java stdout."""
    total_runtime_match = re.search(r'Total runtime:\s*([\d.]+)\s*ms', stdout)
    if total_runtime_match is None:
        total_runtime_match = re.search(r'Completed full run in\s+([\d.]+)\s*ms', stdout)

    timing_window_match = re.search(
        r'Inciso 1\.1 runtime \(first\s+([\d.]+)\s+simulated seconds\):\s*([\d.]+)\s*ms',
        stdout,
    )
    if timing_window_match is None:
        timing_window_match = re.search(
            r'Completed first\s+([\d.]+)\s+simulated seconds in\s+([\d.]+)\s*ms',
            stdout,
        )

    return {
        "total_runtime_ms": float(total_runtime_match.group(1)) if total_runtime_match else None,
        "measurement_window_s": float(timing_window_match.group(1)) if timing_window_match else None,
        "measurement_window_ms": float(timing_window_match.group(2)) if timing_window_match else None,
    }


def run_batch(
    n_values,
    runs_per_n,
    seed_base=42,
    t_final=DEFAULT_TIMING_WINDOW_1_1,
    timing_filename="timing.txt",
    use_measurement_window=False,
):
    """
    Run the Java simulation for various N values and collect timing data.
    Writes timing.txt to data/ directory.
    """
    os.makedirs("data", exist_ok=True)
    timing_file = os.path.join("data", timing_filename)
    measurement_window = t_final if use_measurement_window else min(DEFAULT_TIMING_WINDOW_1_1, t_final)
    stored_metric = "measurement_window" if use_measurement_window else "total_runtime"
    write_output = not use_measurement_window

    with open(timing_file, 'w') as tf:
        tf.write("# N time_ms\n")
        tf.write(
            f"# runs_per_n={runs_per_n}, seed_base={seed_base}, "
            f"sim_t_final={t_final}, measurement_window_s={measurement_window}, "
            f"stored_metric={stored_metric}, write_output={'true' if write_output else 'false'}\n"
        )
    
    for N in n_values:
        for run in range(runs_per_n):
            seed = seed_base + run
            print(f"Running N={N}, run {run+1}/{runs_per_n} (seed={seed})...")
            
            cmd = [
                "java", "-cp", "engine/target/classes",
                "ar.edu.itba.sds.tp3.EventDrivenSimulation",
                "-N", str(N), "-seed", str(seed), "-t_final", str(t_final)
            ]
            if use_measurement_window:
                cmd.extend(["--timing-window", str(measurement_window)])
                cmd.append("--no-output")
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=".")
            
            if result.returncode != 0:
                print(f"  ERROR: {result.stderr}")
                continue
            
            timings = parse_runtime_report(result.stdout)
            if use_measurement_window:
                measured_window_s = timings["measurement_window_s"]
                time_ms = timings["measurement_window_ms"]
            else:
                measured_window_s = None
                time_ms = timings["total_runtime_ms"]
            
            if time_ms is not None:
                with open(timing_file, 'a') as tf:
                    tf.write(f"{N} {time_ms:.2f}\n")
                if use_measurement_window and measured_window_s is not None:
                    print(f"  Completed: {time_ms:.2f} ms hasta t={measured_window_s:g} s simulados")
                else:
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
        return list(range(100, args.n_max + 1, 50))

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
                        help='Use the inciso 1.1 sweep and store timings up to --t-final simulated seconds')
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
    parser.add_argument('--t-final', type=float, default=DEFAULT_TIMING_WINDOW_1_1,
                        help='Simulation final time t_f in seconds')
    args = parser.parse_args()

    t_final = args.t_final
    timing_filename = "timing_1_1.txt" if args.for_1_1 else "timing.txt"
    run_batch(
        build_n_values(args),
        args.runs,
        args.seed,
        t_final,
        timing_filename,
        use_measurement_window=args.for_1_1,
    )


if __name__ == '__main__':
    main()
