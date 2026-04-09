"""
animate_system1.py - Animation for System 1: Scanning rate simulation
Visualizes the particles moving inside the circular enclosure with the central obstacle.
FRESH particles are shown in green, USED particles in violet.

Usage:
    python graphics/animate_system1.py data/sim_100N_20260410_1530.txt
    python graphics/animate_system1.py data/sim_100N_20260410_1530.txt --save output.gif
    python graphics/animate_system1.py data/sim_100N_20260410_1530.txt --interval 50
"""

import sys
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Circle


def parse_simulation_file(filepath):
    """Parse simulation file and return metadata + snapshots."""
    metadata = {}
    snapshots = []
    
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
                particles.append({
                    'id': int(parts[0]),
                    'x': float(parts[1]),
                    'y': float(parts[2]),
                    'vx': float(parts[3]),
                    'vy': float(parts[4]),
                    'state': parts[5]
                })
            snapshots.append((time, particles))
        i += 1
    
    return metadata, snapshots


def animate(filepath, interval_ms=50, save_path=None, skip_frames=1):
    """Create animation of the simulation."""
    metadata, snapshots = parse_simulation_file(filepath)
    
    N = int(metadata['N'])
    R_enclosure = float(metadata['R_enclosure'])
    r0 = float(metadata['r0'])
    r_particle = float(metadata['r'])
    
    # Subsample frames for smooth animation
    snapshots = snapshots[::skip_frames]
    n_frames = len(snapshots)
    
    print(f"Animating {n_frames} frames, N={N}")
    
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_xlim(-R_enclosure - 2, R_enclosure + 2)
    ax.set_ylim(-R_enclosure - 2, R_enclosure + 2)
    ax.set_aspect('equal')
    ax.set_facecolor('#1a1a2e')
    fig.patch.set_facecolor('#0f0f23')
    
    # Draw enclosure
    enclosure = Circle((0, 0), R_enclosure, fill=False,
                        edgecolor='#e0e0e0', linewidth=2, linestyle='-')
    ax.add_patch(enclosure)
    
    # Draw obstacle
    obstacle = Circle((0, 0), r0, fill=True,
                       facecolor='#555555', edgecolor='#888888', linewidth=1.5)
    ax.add_patch(obstacle)
    
    # Initialize scatter plots for fresh and used particles
    scatter_fresh = ax.scatter([], [], s=15, c='#00E676', edgecolors='#00C853',
                                linewidths=0.5, zorder=3, label='Frescas')
    scatter_used = ax.scatter([], [], s=15, c='#D500F9', edgecolors='#AA00FF',
                               linewidths=0.5, zorder=3, label='Usadas')
    
    time_text = ax.text(0.02, 0.98, '', transform=ax.transAxes,
                        color='white', fontsize=14, va='top',
                        bbox=dict(facecolor='black', alpha=0.7, edgecolor='none'))
    
    count_text = ax.text(0.02, 0.90, '', transform=ax.transAxes,
                         color='white', fontsize=11, va='top',
                         bbox=dict(facecolor='black', alpha=0.7, edgecolor='none'))
    
    ax.set_title('Sistema 1: Scanning Rate', color='white', fontsize=16, fontweight='bold')
    ax.set_xlabel('x [m]', color='white', fontsize=12)
    ax.set_ylabel('y [m]', color='white', fontsize=12)
    ax.tick_params(colors='white')
    ax.legend(loc='upper right', fontsize=10, facecolor='#2a2a4a',
              edgecolor='white', labelcolor='white')
    
    def update(frame_idx):
        t, particles = snapshots[frame_idx]
        
        fresh_x, fresh_y = [], []
        used_x, used_y = [], []
        
        for p in particles:
            if p['state'] == 'F':
                fresh_x.append(p['x'])
                fresh_y.append(p['y'])
            else:
                used_x.append(p['x'])
                used_y.append(p['y'])
        
        scatter_fresh.set_offsets(np.c_[fresh_x, fresh_y] if fresh_x else np.empty((0, 2)))
        scatter_used.set_offsets(np.c_[used_x, used_y] if used_x else np.empty((0, 2)))
        
        time_text.set_text(f't = {t:.4f} s')
        count_text.set_text(f'Frescas: {len(fresh_x)} | Usadas: {len(used_x)}')
        
        return scatter_fresh, scatter_used, time_text, count_text
    
    anim = FuncAnimation(fig, update, frames=n_frames,
                         interval=interval_ms, blit=False, repeat=True)
    
    if save_path:
        suffix = os.path.splitext(save_path)[1].lower()
        if suffix == '.gif':
            fps = 1000 / max(interval_ms, 1)
            anim.save(save_path, writer=PillowWriter(fps=fps))
        elif suffix in {'.mp4', '.m4v'}:
            anim.save(save_path, writer='ffmpeg')
        else:
            anim.save(save_path)
        print(f"Animation saved to: {save_path}")
        plt.close(fig)
    else:
        plt.tight_layout()
        plt.show()


def main():
    parser = argparse.ArgumentParser(description='Animation for System 1')
    parser.add_argument('file', help='Path to simulation output file')
    parser.add_argument('--interval', type=int, default=50,
                        help='Frame interval in ms (default: 50)')
    parser.add_argument('--save', default=None,
                        help='Save animation to file (e.g., output.gif)')
    parser.add_argument('--skip', type=int, default=1,
                        help='Frame skip factor (default: 1)')
    args = parser.parse_args()
    
    animate(args.file, args.interval, args.save, args.skip)


if __name__ == '__main__':
    main()
