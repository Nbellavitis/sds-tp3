"""
plot_metrics.py - Main post-processing orchestrator for System 1
Reads simulation output files from data/ and generates all required plots.

Usage:
    python graphics/plot_metrics.py <path_to_sim_file>      (single file: 1.3 + 1.4)
    python graphics/plot_metrics.py data/                    (directory: 1.1 + 1.2 + 1.3 + 1.4)
"""

import sys
import os
import glob

from plot_execution_time import plot_execution_time
from plot_scanning_rate import parse_simulation_file, plot_scanning_rate
from plot_fraction_used import plot_fraction_used, plot_fest_vs_N
from plot_radial_profiles import plot_radial_profiles, plot_at_S2_vs_N


def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_metrics.py <sim_file_or_directory>")
        print("  Single file: generates plots 1.3 and 1.4")
        print("  Directory:   generates all plots 1.1-1.4")
        sys.exit(1)
    
    path = sys.argv[1]
    output_dir = "graphics/output"
    
    if os.path.isfile(path):
        print(f"Processing single file: {path}")
        metadata, snapshots, events = parse_simulation_file(path)
        N = int(metadata['N'])
        print(f"  N={N}, snapshots={len(snapshots)}, events={len(events)}")
        
        plot_fraction_used(metadata, snapshots, path, output_dir)
        plot_radial_profiles(metadata, snapshots, path, output_dir)
        
    elif os.path.isdir(path):
        print(f"Processing directory: {path}")
        files = sorted(glob.glob(os.path.join(path, 'sim_*.txt')))
        
        if not files:
            print("No simulation files found in directory.")
            sys.exit(1)
        
        print(f"  Found {len(files)} simulation files")
        
        # Group by N
        files_by_N = {}
        for fpath in files:
            meta, snaps, events = parse_simulation_file(fpath)
            n = int(meta['N'])
            if n not in files_by_N:
                files_by_N[n] = []
            files_by_N[n].append((fpath, meta, snaps, events))
        
        print(f"  N values: {sorted(files_by_N.keys())}")
        
        # 1.1 - Execution time
        print("\n=== Inciso 1.1: Execution Time ===")
        plot_execution_time(output_dir=output_dir)
        
        # 1.2 - Scanning rate
        print("\n=== Inciso 1.2: Scanning Rate ===")
        plot_scanning_rate(files_by_N, output_dir)
        
        # 1.3 - Fraction used
        print("\n=== Inciso 1.3: Fraction Used ===")
        for N in sorted(files_by_N.keys()):
            fpath, meta, snaps, _ = files_by_N[N][0]
            plot_fraction_used(meta, snaps, fpath, output_dir)
        if len(files_by_N) > 1:
            plot_fest_vs_N(files_by_N, output_dir)
        
        # 1.4 - Radial profiles
        print("\n=== Inciso 1.4: Radial Profiles ===")
        for N in sorted(files_by_N.keys()):
            fpath, meta, snaps, _ = files_by_N[N][0]
            plot_radial_profiles(meta, snaps, fpath, output_dir)
        if len(files_by_N) > 1:
            plot_at_S2_vs_N(files_by_N, output_dir)
    else:
        print(f"Error: '{path}' is not a valid file or directory.")
        sys.exit(1)
    
    print("\nAll plots generated successfully!")


if __name__ == '__main__':
    main()
