"""
animate_system1_mru.py - Snapshot animation for System 1

Builds a discrete animation directly from saved snapshots in sim_*.txt.
No interpolation is performed between snapshots.

Usage:
    python graphics/animate_system1_mru.py data/sim_100N_20260410_1530.txt
    python graphics/animate_system1_mru.py data/sim_100N_20260410_1530.txt --save graphics/output/system1.gif
    python graphics/animate_system1_mru.py data/sim_100N_20260410_1530.txt --fps 30 --speed 2.0
"""

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Circle


FRESH_FACE = "#28C76F"
FRESH_EDGE = "#1A8F4B"
USED_FACE = "#D633FF"
USED_EDGE = "#8F1FB3"


def parse_simulation_file(filepath):
    """Parse a System 1 output file into dense numpy arrays."""
    metadata = {}
    snapshot_times = []
    xs = []
    ys = []
    vxs = []
    vys = []
    used_flags = []

    with open(filepath, "r") as f:
        lines = f.readlines()

    i = 0
    total = len(lines)

    while i < total and lines[i].startswith("#"):
        line = lines[i].strip().lstrip("# ")
        if "N=" in line:
            for part in line.split():
                if "=" not in part:
                    continue
                key, value = part.split("=", 1)
                try:
                    metadata[key] = float(value)
                except ValueError:
                    metadata[key] = value
        i += 1

    N = int(metadata.get("N", 0))

    while i < total:
        line = lines[i].strip()
        if not line:
            i += 1
            continue

        if line.startswith("S "):
            snapshot_times.append(float(line.split()[1]))
            x = np.empty(N, dtype=float)
            y = np.empty(N, dtype=float)
            vx = np.empty(N, dtype=float)
            vy = np.empty(N, dtype=float)
            used = np.zeros(N, dtype=bool)

            for _ in range(N):
                i += 1
                if i >= total:
                    break
                parts = lines[i].strip().split()
                if len(parts) < 6:
                    continue
                pid = int(parts[0])
                x[pid] = float(parts[1])
                y[pid] = float(parts[2])
                vx[pid] = float(parts[3])
                vy[pid] = float(parts[4])
                used[pid] = parts[5] == "U"

            xs.append(x)
            ys.append(y)
            vxs.append(vx)
            vys.append(vy)
            used_flags.append(used)

        i += 1

    if not snapshot_times:
        raise ValueError(f"No snapshots found in {filepath}")

    return {
        "metadata": metadata,
        "times": np.array(snapshot_times, dtype=float),
        "x": np.stack(xs),
        "y": np.stack(ys),
        "vx": np.stack(vxs),
        "vy": np.stack(vys),
        "used": np.stack(used_flags),
    }


def select_snapshot_indices(snapshot_count, speed):
    """Select discrete snapshot indices, optionally skipping with speed>1."""
    stride = max(1, int(round(speed)))
    return np.arange(0, snapshot_count, stride, dtype=int), stride


def animate(filepath, fps=30, speed=1.0, save_path=None, dpi=120):
    """Create and optionally save a snapshot-only animation."""
    if fps <= 0:
        raise ValueError("fps must be positive")
    if speed <= 0:
        raise ValueError("speed must be positive")

    data = parse_simulation_file(filepath)
    metadata = data["metadata"]

    N = int(metadata["N"])
    R_enclosure = float(metadata["R_enclosure"])
    r0 = float(metadata["r0"])
    r_particle = float(metadata["r"])
    snapshot_every_events = int(metadata.get("snapshot_every_events", 1))
    frame_indices, snapshot_stride = select_snapshot_indices(len(data["times"]), speed)
    frame_times = data["times"][frame_indices]
    n_frames = len(frame_indices)
    playback_seconds = n_frames / fps
    print(f"Animating {os.path.basename(filepath)}")
    print(f"  N={N}, snapshots={len(data['times'])}, frames={n_frames}")
    if snapshot_every_events > 1:
        print(f"  Note: file was saved every {snapshot_every_events} collisions (discrete snapshot playback).")
    if snapshot_stride > 1:
        print(f"  Using one every {snapshot_stride} snapshots due to speed={speed:.2f}")
    print(f"  simulated time: {frame_times[0]:.3f} s -> {frame_times[-1]:.3f} s")
    print(f"  playback time: ~{playback_seconds:.2f} s at {fps} fps")

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_aspect("equal")
    ax.set_xlim(-R_enclosure - 2.0, R_enclosure + 2.0)
    ax.set_ylim(-R_enclosure - 2.0, R_enclosure + 2.0)
    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")

    enclosure = Circle((0.0, 0.0), R_enclosure, fill=False, edgecolor="black", linewidth=1.8)
    obstacle = Circle((0.0, 0.0), r0, facecolor="black", edgecolor="black", linewidth=1.0)
    ax.add_patch(enclosure)
    ax.add_patch(obstacle)

    particles = []
    for _ in range(N):
        particle = Circle((0.0, 0.0), r_particle, linewidth=0.6, zorder=3)
        ax.add_patch(particle)
        particles.append(particle)

    # Formato limpio tipo consigna: sin ejes ni etiquetas superpuestas.
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    interval_ms = max(1, int(round(1000.0 / fps)))

    def update(frame_idx):
        snapshot_idx = frame_indices[frame_idx]
        x = data["x"][snapshot_idx]
        y = data["y"][snapshot_idx]
        used = data["used"][snapshot_idx]

        for pid, particle in enumerate(particles):
            particle.center = (x[pid], y[pid])
            if used[pid]:
                particle.set_facecolor(USED_FACE)
                particle.set_edgecolor(USED_EDGE)
            else:
                particle.set_facecolor(FRESH_FACE)
                particle.set_edgecolor(FRESH_EDGE)

        return particles

    anim = FuncAnimation(fig, update, frames=n_frames, interval=interval_ms, blit=False, repeat=True)

    if save_path:
        suffix = os.path.splitext(save_path)[1].lower()
        os.makedirs(os.path.dirname(save_path) or ".", exist_ok=True)
        if suffix == ".gif":
            anim.save(save_path, writer=PillowWriter(fps=fps), dpi=dpi)
        elif suffix in {".mp4", ".m4v"}:
            anim.save(save_path, writer="ffmpeg", dpi=dpi, fps=fps)
        else:
            anim.save(save_path, dpi=dpi)
        print(f"Animation saved to: {save_path}")
        plt.close(fig)
    else:
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Animate System 1 using saved snapshots only (no interpolation)."
    )
    parser.add_argument("file", help="Path to one simulation output file")
    parser.add_argument("--fps", type=int, default=30, help="Playback frames per second")
    parser.add_argument(
        "--speed",
        type=float,
        default=1.0,
        help="Snapshot skip control (1=all snapshots, 2~one every 2 snapshots, etc.)",
    )
    parser.add_argument("--save", default=None, help="Optional output path, e.g. graphics/output/system1.gif")
    parser.add_argument("--dpi", type=int, default=120, help="Output DPI when saving")
    args = parser.parse_args()

    animate(args.file, fps=args.fps, speed=args.speed, save_path=args.save, dpi=args.dpi)


if __name__ == "__main__":
    main()
