#!/usr/bin/env python3
"""Prepare and run the UTI89 decision-tree stage."""

import argparse
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from uti89_pipeline.config import load_config  # noqa: E402
from uti89_pipeline.decision_tree import (  # noqa: E402
    format_command,
    prepare_decision_tree_run,
    run_decision_tree,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="Pipeline YAML config.")
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Prepare the decision-tree run directory but do not write decisions.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Write a decision-tree plan without running external searches or writing decisions.",
    )
    parser.add_argument(
        "--unlock",
        action="store_true",
        help="Accepted for compatibility with older configs; this is a no-op.",
    )
    args = parser.parse_args()

    config = load_config(args.config)
    run = prepare_decision_tree_run(config)

    print(f"Prepared decision-tree run directory: {run.work_dir}")
    print(f"Intermediate directory: {run.pipeline_files}")
    print(f"Expected decisions file: {run.decisions}")
    print("Decision-tree command:")
    print(format_command(run.command, run.work_dir))

    if args.prepare_only:
        return 0
    return run_decision_tree(run, dry_run=args.dry_run, unlock=args.unlock)


if __name__ == "__main__":
    raise SystemExit(main())
