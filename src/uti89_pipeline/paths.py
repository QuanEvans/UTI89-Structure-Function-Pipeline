"""Path helpers for pipeline stages."""

from pathlib import Path
from typing import Any, Dict


def run_dir(config: Dict[str, Any]) -> Path:
    return Path(config["run"]["work_dir"]).expanduser().resolve()


def decision_file(config: Dict[str, Any]) -> Path:
    root = run_dir(config)
    structure_cfg = config.get("structure_prediction", {})
    return Path(
        structure_cfg.get(
            "decision_file",
            str(root / "decision_tree_intermediates" / "decide" / "decisions.txt"),
        )
    ).expanduser().resolve()


def structure_output_dir(config: Dict[str, Any]) -> Path:
    root = run_dir(config)
    structure_cfg = config.get("structure_prediction", {})
    return Path(structure_cfg.get("output_dir", str(root / "structures"))).expanduser().resolve()


def function_output_dir(config: Dict[str, Any]) -> Path:
    root = run_dir(config)
    function_cfg = config.get("function_prediction", {})
    return Path(function_cfg.get("output_dir", str(root / "functions"))).expanduser().resolve()

