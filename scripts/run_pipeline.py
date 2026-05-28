#!/usr/bin/env python3
"""Run the UTI89 pipeline stages in order from one config file."""

import argparse
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from uti89_pipeline.config import load_config  # noqa: E402
from uti89_pipeline.decision_tree import prepare_decision_tree_run, run_decision_tree  # noqa: E402
from uti89_pipeline.decisions import read_decision_map  # noqa: E402
from uti89_pipeline.function import prepare_function_predictions  # noqa: E402
from uti89_pipeline.paths import decision_file, function_output_dir, structure_output_dir  # noqa: E402
from uti89_pipeline.structure import prepare_structure_predictions  # noqa: E402
from uti89_pipeline.wait import wait_for_files  # noqa: E402


STAGES = ["decision", "structure", "function"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="Pipeline YAML config.")
    parser.add_argument(
        "--work-dir",
        help="Override run.work_dir from the config. Useful for repeatable smoke tests.",
    )
    parser.add_argument(
        "--submit",
        action="store_true",
        help="Run generated structure/function jobs with the configured execution backend.",
    )
    parser.add_argument(
        "--wait",
        action="store_true",
        help="Wait for structure/function outputs after submitting jobs.",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=int,
        default=86400,
        help="Maximum wait time per stage when --wait is used.",
    )
    parser.add_argument(
        "--poll-seconds",
        type=int,
        default=60,
        help="Polling interval when --wait is used.",
    )
    parser.add_argument(
        "--start-at",
        choices=STAGES,
        default="decision",
        help="First stage to run.",
    )
    parser.add_argument(
        "--stop-after",
        choices=STAGES,
        default="function",
        help="Last stage to run.",
    )
    parser.add_argument(
        "--decision-dry-run",
        action="store_true",
        help="Dry-run only the native decision-tree stage.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Regenerate structure/function jobs even if outputs already exist.",
    )
    args = parser.parse_args()

    config = load_config(args.config)
    if args.work_dir:
        _apply_work_dir_override(config, args.work_dir)
    selected = _selected_stages(args.start_at, args.stop_after)

    if "decision" in selected:
        if _uses_user_models(config):
            print("== Decision tree ==")
            print("Skipping decision tree because structure_prediction.user_models is configured.")
        else:
            print("== Decision tree ==")
            run = prepare_decision_tree_run(config)
            rc = run_decision_tree(run, dry_run=args.decision_dry_run)
            if rc != 0:
                return rc
            if args.decision_dry_run:
                return 0
            _require_file(run.decisions, "decision-tree output")
            print("Decision output: {}".format(run.decisions))

    if "structure" in selected:
        print("== Structure prediction ==")
        result = prepare_structure_predictions(config, submit=args.submit, force=args.force)
        print("Structure output directory: {}".format(result.output_dir))
        print("Prepared {} structure actions.".format(len(result.prepared)))
        _print_skips(result.skipped)
        if args.submit and args.wait:
            missing = wait_for_files(
                _structure_expected_outputs(config),
                timeout_seconds=args.timeout_seconds,
                poll_seconds=args.poll_seconds,
            )
            if missing:
                _print_missing("Timed out waiting for structure outputs", missing)
                return 1

    if "function" in selected:
        print("== Function prediction ==")
        result = prepare_function_predictions(config, submit=args.submit, force=args.force)
        print("Function output directory: {}".format(result.output_dir))
        print("Prepared {} function jobs.".format(len(result.prepared)))
        _print_skips(result.skipped)
        if args.submit and args.wait:
            missing = wait_for_files(
                _function_expected_outputs(config),
                timeout_seconds=args.timeout_seconds,
                poll_seconds=args.poll_seconds,
            )
            if missing:
                _print_missing("Timed out waiting for function outputs", missing)
                return 1

    return 0


def _selected_stages(start_at: str, stop_after: str):
    start = STAGES.index(start_at)
    stop = STAGES.index(stop_after)
    if start > stop:
        raise SystemExit("--start-at must be before or equal to --stop-after")
    return STAGES[start : stop + 1]


def _apply_work_dir_override(config, work_dir: str) -> None:
    """Point derived stage paths at the overridden run directory."""
    config.setdefault("run", {})["work_dir"] = work_dir
    for section, keys in {
        "structure_prediction": ("output_dir", "decision_file"),
        "function_prediction": ("structure_dir", "output_dir", "decision_file"),
    }.items():
        section_cfg = config.get(section)
        if isinstance(section_cfg, dict):
            for key in keys:
                section_cfg.pop(key, None)


def _uses_user_models(config) -> bool:
    user_models = config.get("structure_prediction", {}).get("user_models", {})
    return isinstance(user_models, dict) and bool(user_models) and user_models.get("enabled", True)


def _structure_expected_outputs(config):
    root = structure_output_dir(config)
    return [root / protein_id / "model1.pdb" for protein_id in _structure_ids(config)]


def _function_expected_outputs(config):
    root = function_output_dir(config)
    return [root / protein_id / "consensus.tsv" for protein_id in _function_input_ids(config)]


def _structure_ids(config):
    decisions = decision_file(config)
    if decisions.is_file():
        return sorted(read_decision_map(decisions).keys())
    root = structure_output_dir(config)
    return _ids_from_structure_dir(root)


def _function_input_ids(config):
    decisions = decision_file(config)
    if decisions.is_file():
        return sorted(read_decision_map(decisions).keys())
    root = _function_structure_dir(config)
    return _ids_from_structure_dir(root)


def _function_structure_dir(config):
    run_root = Path(config["run"]["work_dir"]).expanduser().resolve()
    structure_cfg = config.get("structure_prediction", {})
    function_cfg = config.get("function_prediction", {})
    return Path(
        function_cfg.get(
            "structure_dir",
            structure_cfg.get("output_dir", str(run_root / "structures")),
        )
    ).expanduser().resolve()


def _ids_from_structure_dir(root):
    if not root.is_dir():
        return []
    return sorted(
        path.name
        for path in root.iterdir()
        if path.is_dir() and (path / "model1.pdb").exists() and (path / "seq.fasta").exists()
    )


def _require_file(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size < 1:
        raise FileNotFoundError("Missing {}: {}".format(label, path))


def _print_skips(skipped) -> None:
    if not skipped:
        return
    print("Skipped {} entries.".format(len(skipped)))
    for entry in skipped[:20]:
        print(entry)
    if len(skipped) > 20:
        print("... {} more".format(len(skipped) - 20))


def _print_missing(message: str, missing) -> None:
    print(message)
    for path in missing[:20]:
        print(path)
    if len(missing) > 20:
        print("... {} more".format(len(missing) - 20))


if __name__ == "__main__":
    raise SystemExit(main())
