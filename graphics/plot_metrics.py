"""
plot_metrics.py - Main post-processing script for System 1: Scanning Rate
Reads the simulation output file from data/ and generates all required plots.

Usage:
    python graphics/plot_metrics.py <path_to_sim_file>
    python graphics/plot_metrics.py data/sim_100N_20260410_1530.txt

Generates: All plots for incisos 1.1 through 1.4
"""

import sys
import os
import numpy as np
from pathlib import Path

# Import sub-modules
from plot_execution_time import plot_execution_time
from plot_scanning_rate import plot_scanning_rate
from plot_fraction_used import plot_fraction_used
from plot_radial_profiles import plot_radial_profiles


def parse_simulation_file(filepath):
    """
    Parse the simulation output file.
    
    Returns:
        metadata: dict with N, L, R_enclosure, r0, r, m, v0, t_final
        snapshots: list of (time, particles_data) where particles_data is
                   a list of dicts with id, x, y, vx, vy, state
    """
    metadata = {}
    snapshots = []
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    i = 0
    total = len(lines)
    
    # Parse header/metadata
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
    
    # Parse snapshots
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
                particles.append({
                    'id': int(parts[0]),
                    'x': float(parts[1]),
                    'y': float(parts[2]),
                    'vx': float(parts[3]),
                    'vy': float(parts[4]),
                    'state': parts[5]  # 'F' or 'U'
                })
            snapshots.append((time, particles))
        i += 1
    
    return metadata, snapshots


def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_metrics.py <sim_file_or_directory>")
        print("  Single file: generates plots 1.3 and 1.4")
        print("  Directory with multiple files: generates all plots 1.1-1.4")
        sys.exit(1)
    
    path = sys.argv[1]
    
    if os.path.isfile(path):
        print(f"Processing single file: {path}")
        metadata, snapshots = parse_simulation_file(path)
        print(f"  N={int(metadata['N'])}, snapshots={len(snapshots)}")
        
        # Generate individual plots
        plot_fraction_used(metadata, snapshots, path)
        plot_radial_profiles(metadata, snapshots, path)
        
    elif os.path.isdir(path):
        print(f"Processing directory: {path}")
        files = sorted([os.path.join(path, f) for f in os.listdir(path)
                        if f.startswith('sim_') and f.endswith('.txt')])
        
        if not files:
            print("No simulation files found in directory.")
            sys.exit(1)
        
        print(f"  Found {len(files)} simulation files")
        
        # Group files by N value
        files_by_N = {}
        for fpath in files:
            meta, snaps = parse_simulation_file(fpath)
            n = int(meta['N'])
            if n not in files_by_N:
                files_by_N[n] = []
            files_by_N[n].append((fpath, meta, snaps))
        
        # Generate all plots
        plot_execution_time(files_by_N)
        plot_scanning_rate(files_by_N)
        
        # For fraction used and radial profiles, use the largest N
        max_N = max(files_by_N.keys())
        fpath, meta, snaps = files_by_N[max_N][0]
        plot_fraction_used(meta, snaps, fpath)
        plot_radial_profiles(meta, snaps, fpath)
    else:
        print(f"Error: '{path}' is not a valid file or directory.")
        sys.exit(1)
    
    print("\nAll plots generated successfully!")


if __name__ == '__main__':
    main()
