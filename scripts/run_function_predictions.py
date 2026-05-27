#!/usr/bin/env python3
"""Prepare and optionally submit function-prediction jobs."""

import argparse
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from uti89_pipeline.config import load_config  # noqa: E402
from uti89_pipeline.function import prepare_function_predictions  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="Pipeline YAML config.")
    parser.add_argument(
        "--submit",
        action="store_true",
        help="Run generated scripts with the configured execution backend.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Prepare jobs even when consensus.tsv already exists.",
    )
    args = parser.parse_args()

    config = load_config(args.config)
    result = prepare_function_predictions(config, submit=args.submit, force=args.force)
    print("Function output directory: {}".format(result.output_dir))
    print("Prepared {} function jobs.".format(len(result.prepared)))
    for path in result.prepared[:20]:
        print(path)
    if len(result.prepared) > 20:
        print("... {} more".format(len(result.prepared) - 20))
    if result.skipped:
        print("Skipped {} entries.".format(len(result.skipped)))
        for entry in result.skipped[:20]:
            print(entry)
        if len(result.skipped) > 20:
            print("... {} more".format(len(result.skipped) - 20))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
