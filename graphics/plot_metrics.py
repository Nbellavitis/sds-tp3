"""
plot_metrics.py - Main post-processing orchestrator for System 1
Reads simulation output files from data/ and generates all required plots.

Usage:
    python graphics/plot_metrics.py <path_to_sim_file>      (single file: 1.3 + 1.4)
    python graphics/plot_metrics.py data/                    (directory: 1.1 + 1.2 + 1.3 + 1.4)
"""

import sys
import os

from analysis_cache import group_entries_by_N, load_analysis_entries, load_analysis_file
from plot_execution_time import plot_execution_time
from plot_scanning_rate import plot_scanning_rate
from plot_fraction_used import (
    plot_fraction_used,
    plot_fraction_used_realizations,
    plot_t_est_vs_N,
    plot_f_est_vs_N,
)
from plot_radial_profiles import plot_radial_profiles, plot_radial_profiles_ensemble, plot_at_S2_vs_N


DISPLAY_PROFILE_N_VALUES = {100, 300, 500, 700}


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
        entry = load_analysis_file(path)
        metadata = entry["metadata"]
        N = int(metadata['N'])
        print(f"  N={N}, F_u samples={len(entry['fu']['times'])}, C_fc samples={len(entry['cfc']['times'])}")
        
        plot_fraction_used(entry=entry, output_dir=output_dir)
        plot_radial_profiles(entry=entry, output_dir=output_dir)
        
    elif os.path.isdir(path):
        print(f"Processing directory: {path}")
        entries = load_analysis_entries(path)

        if not entries:
            print("No simulation files found in directory.")
            sys.exit(1)
        
        print(f"  Found {len(entries)} simulation files")
        files_by_N = group_entries_by_N(entries)
        selected_profile_ns = [N for N in sorted(files_by_N.keys()) if N in DISPLAY_PROFILE_N_VALUES]
        selected_files_by_N = {N: files_by_N[N] for N in selected_profile_ns}

        print(f"  N values: {sorted(files_by_N.keys())}")
        print(f"  N mostrados en figuras por realizacion/perfil: {selected_profile_ns}")
        
        # 1.1 - Execution time
        print("\n=== Inciso 1.1: Execution Time ===")
        plot_execution_time(data_dir=path, output_dir=output_dir)
        
        # 1.2 - Scanning rate
        print("\n=== Inciso 1.2: Scanning Rate ===")
        plot_scanning_rate(files_by_N, output_dir)
        
        # 1.3 - Fraction used
        print("\n=== Inciso 1.3: Fraction Used ===")
        for N in selected_profile_ns:
            plot_fraction_used_realizations(files_by_N[N], output_dir)
        plot_t_est_vs_N(selected_files_by_N, output_dir)
        plot_f_est_vs_N(selected_files_by_N, output_dir)

        # 1.4 - Radial profiles
        print("\n=== Inciso 1.4: Radial Profiles ===")
        for N in selected_profile_ns:
            plot_radial_profiles_ensemble(files_by_N[N], output_dir)
        if len(selected_files_by_N) > 1:
            plot_at_S2_vs_N(selected_files_by_N, output_dir)
    else:
        print(f"Error: '{path}' is not a valid file or directory.")
        sys.exit(1)
    
    print("\nAll plots generated successfully!")


if __name__ == '__main__':
    main()
