"""
analysis_cache.py - Lightweight Python loader for Java-built analysis caches.

The heavy parsing and per-file metric extraction are delegated to the
SimulationAnalysisCacheBuilder Java class. Python only ensures caches exist,
loads the compact JSON payload, and plots from there.
"""

import gzip
import json
import os
import subprocess
from pathlib import Path


CACHE_SUFFIX = ".analysis.v2.json.gz"
ANALYZER_CLASS = "ar.edu.itba.sds.tp3.SimulationAnalysisCacheBuilder"
REPO_ROOT = Path(__file__).resolve().parent.parent
ENGINE_CLASSES = REPO_ROOT / "engine" / "target" / "classes"
ANALYZER_CLASS_FILE = ENGINE_CLASSES / "ar" / "edu" / "itba" / "sds" / "tp3" / "SimulationAnalysisCacheBuilder.class"
ANALYZER_SOURCE_FILE = REPO_ROOT / "engine" / "src" / "main" / "java" / "ar" / "edu" / "itba" / "sds" / "tp3" / "SimulationAnalysisCacheBuilder.java"
DEFAULT_CACHE_DIR = REPO_ROOT / "data" / "cache"


def list_sim_files(path):
    """Return sorted sim_*.txt files from a file or directory path."""
    path = Path(path)
    if path.is_file():
        return [path.resolve()]

    files = sorted(
        p.resolve()
        for p in path.iterdir()
        if p.is_file() and p.name.startswith("sim_") and p.suffix == ".txt"
    )
    return files


def cache_path_for(sim_path, cache_dir=DEFAULT_CACHE_DIR):
    """Map one simulation file path to its cache path."""
    sim_path = Path(sim_path)
    cache_dir = Path(cache_dir)
    return cache_dir / f"{sim_path.stem}{CACHE_SUFFIX}"


def backend_needs_compile():
    """Check whether the Java analyzer class is missing or stale."""
    if not ANALYZER_CLASS_FILE.exists():
        return True
    if ANALYZER_SOURCE_FILE.exists() and ANALYZER_SOURCE_FILE.stat().st_mtime > ANALYZER_CLASS_FILE.stat().st_mtime:
        return True
    return False


def ensure_backend_compiled():
    """Compile the Java backend if needed."""
    if not backend_needs_compile():
        return

    subprocess.run(
        ["mvn", "-f", "engine/pom.xml", "compile"],
        cwd=REPO_ROOT,
        check=True,
    )


def cache_is_fresh(sim_path, cache_path):
    """Check whether a cache file exists and is newer than its source file."""
    sim_path = Path(sim_path)
    cache_path = Path(cache_path)
    if not cache_path.exists():
        return False
    return cache_path.stat().st_mtime >= sim_path.stat().st_mtime


def ensure_analysis_cache(path, cache_dir=DEFAULT_CACHE_DIR):
    """
    Ensure the cache exists for one simulation file or for all files in a directory.
    """
    path = Path(path).resolve()
    cache_dir = Path(cache_dir).resolve()
    sim_files = list_sim_files(path)

    if not sim_files:
        raise FileNotFoundError(f"No simulation files found in {path}")

    stale_files = [sim for sim in sim_files if not cache_is_fresh(sim, cache_path_for(sim, cache_dir))]
    if not stale_files:
        return

    ensure_backend_compiled()
    cache_dir.mkdir(parents=True, exist_ok=True)

    target = path if path.is_dir() else stale_files[0]
    subprocess.run(
        [
            "java",
            "-cp",
            str(ENGINE_CLASSES),
            ANALYZER_CLASS,
            str(target),
            "--cache-dir",
            str(cache_dir),
        ],
        cwd=REPO_ROOT,
        check=True,
    )


def load_analysis_file(sim_path, cache_dir=DEFAULT_CACHE_DIR):
    """Load one cached analysis entry, building the cache if needed."""
    sim_path = Path(sim_path).resolve()
    ensure_analysis_cache(sim_path, cache_dir)
    cache_path = cache_path_for(sim_path, cache_dir)
    with gzip.open(cache_path, "rt", encoding="utf-8") as f:
        data = json.load(f)
    data["source_path"] = str(sim_path)
    return data


def load_analysis_entries(path, cache_dir=DEFAULT_CACHE_DIR):
    """Load all cached analysis entries for a file or directory."""
    path = Path(path).resolve()
    ensure_analysis_cache(path, cache_dir)
    entries = []
    for sim_path in list_sim_files(path):
        entries.append(load_analysis_file(sim_path, cache_dir))
    return entries


def group_entries_by_N(entries):
    """Group cached entries by particle count N."""
    grouped = {}
    for entry in entries:
        n = int(entry["metadata"]["N"])
        grouped.setdefault(n, []).append(entry)
    return grouped
