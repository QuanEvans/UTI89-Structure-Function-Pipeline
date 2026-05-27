"""Polling helpers for Slurm-produced outputs."""

from pathlib import Path
from typing import Iterable, List
import time


def nonempty_file(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 1


def wait_for_files(
    paths: Iterable[Path],
    timeout_seconds: int,
    poll_seconds: int,
    require_nonempty: bool = True,
) -> List[Path]:
    """Wait for expected files and return any files still missing."""
    expected = list(paths)
    deadline = time.time() + timeout_seconds
    while True:
        missing = [
            path
            for path in expected
            if not (nonempty_file(path) if require_nonempty else path.exists())
        ]
        if not missing:
            return []
        if time.time() >= deadline:
            return missing
        time.sleep(poll_seconds)

