"""
run_animation_sim.py - Runner dedicado para generar salidas de simulacion orientadas a animacion.

No reemplaza run_batch.py ni el flujo de graficos ya existente.
Solo ejecuta el motor Java con snapshot configurable para no perder eventos
entre snapshots cuando se quiera animacion mas fiel.

Uso:
    python graphics/run_animation_sim.py --n 300 --seed 42 --t-final 8.0 --snapshot-every-events 1
    python graphics/run_animation_sim.py --n 300 --runs 3 --seed 42 --snapshot-every-events 2
"""

import argparse
import re
import subprocess
import sys


def extract_output_path(stdout):
    """Extract generated output path from Java stdout, if present."""
    match = re.search(r"Output:\s*(.+)", stdout)
    if not match:
        return None
    value = match.group(1).strip()
    if value.endswith("(--no-output)"):
        return None
    return value


def run_animation_sim(n, seed_base, runs, t_final, snapshot_every_events):
    """Run Java simulation(s) with configurable snapshot frequency for animation."""
    generated_files = []

    for run in range(runs):
        seed = seed_base + run
        cmd = [
            "java",
            "-cp",
            "engine/target/classes",
            "ar.edu.itba.sds.tp3.EventDrivenSimulation",
            "-N",
            str(n),
            "-seed",
            str(seed),
            "-runs",
            "1",
            "-t_final",
            str(t_final),
            "-snapshot_every_events",
            str(snapshot_every_events),
        ]

        print(f"[animation-runner] Running N={n}, seed={seed}, snapshot_every_events={snapshot_every_events}")
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=".")

        if result.returncode != 0:
            print("[animation-runner] Java execution failed:", file=sys.stderr)
            if result.stdout:
                print(result.stdout, file=sys.stderr)
            if result.stderr:
                print(result.stderr, file=sys.stderr)
            raise RuntimeError("Java simulation failed")

        if result.stdout:
            print(result.stdout.strip())

        output_path = extract_output_path(result.stdout)
        if output_path:
            generated_files.append(output_path)

    return generated_files


def main():
    parser = argparse.ArgumentParser(
        description="Run simulation with configurable snapshot frequency for animation files"
    )
    parser.add_argument("--n", type=int, required=True, help="Number of particles")
    parser.add_argument("--seed", type=int, default=42, help="Base seed (run index is added)")
    parser.add_argument("--runs", type=int, default=1, help="Number of runs")
    parser.add_argument("--t-final", type=float, default=5.0, help="Simulation final time")
    parser.add_argument(
        "--snapshot-every-events",
        type=int,
        default=1,
        help="Save one snapshot every K collisions (K>=1). Use 1 for max detail.",
    )
    args = parser.parse_args()

    if args.n <= 0:
        raise ValueError("--n must be positive")
    if args.runs <= 0:
        raise ValueError("--runs must be positive")
    if args.t_final <= 0:
        raise ValueError("--t-final must be positive")
    if args.snapshot_every_events <= 0:
        raise ValueError("--snapshot-every-events must be >= 1")

    files = run_animation_sim(
        n=args.n,
        seed_base=args.seed,
        runs=args.runs,
        t_final=args.t_final,
        snapshot_every_events=args.snapshot_every_events,
    )

    if files:
        print("[animation-runner] Generated files:")
        for fpath in files:
            print(f"  - {fpath}")
    else:
        print("[animation-runner] No output files were detected in stdout.")


if __name__ == "__main__":
    main()

