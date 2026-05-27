"""Prepare function-prediction jobs from structure outputs."""

from pathlib import Path
from typing import Any, Dict, List
import os
import shlex

from .config import get_required
from .decisions import read_decision_map
from .execution import module_load_line, script_suffix, slurm_header, submit_script


class FunctionPreparationResult(object):
    def __init__(self, output_dir: Path, prepared: List[Path], skipped: List[str]) -> None:
        self.output_dir = output_dir
        self.prepared = prepared
        self.skipped = skipped


def prepare_function_predictions(
    config: Dict[str, Any], submit: bool = False, force: bool = False
) -> FunctionPreparationResult:
    """Prepare per-protein StarFunc inputs and Slurm scripts."""
    run_dir = Path(get_required(config, "run.work_dir")).expanduser().resolve()
    function_cfg = config.get("function_prediction", {})
    structure_cfg = config.get("structure_prediction", {})
    structure_dir = Path(
        function_cfg.get(
            "structure_dir",
            structure_cfg.get("output_dir", str(run_dir / "structures")),
        )
    ).expanduser().resolve()
    output_dir = Path(
        function_cfg.get("output_dir", str(run_dir / "functions"))
    ).expanduser().resolve()
    decision_file = Path(
        function_cfg.get(
            "decision_file",
            structure_cfg.get(
                "decision_file",
                str(run_dir / "decision_tree_intermediates" / "decide" / "decisions.txt"),
            ),
        )
    ).expanduser().resolve()
    link_inputs = bool(function_cfg.get("link_inputs", True))

    if not structure_dir.is_dir():
        raise NotADirectoryError("Missing structure directory: {}".format(structure_dir))

    output_dir.mkdir(parents=True, exist_ok=True)
    if decision_file.is_file():
        protein_ids = sorted(read_decision_map(decision_file).keys())
    else:
        protein_ids = _discover_structure_ids(structure_dir)
    prepared = []  # type: List[Path]
    skipped = []  # type: List[str]

    for protein_id in protein_ids:
        src_dir = structure_dir / protein_id
        model = src_dir / "model1.pdb"
        seq = src_dir / "seq.fasta"
        if not _nonempty_file(model):
            skipped.append("{}: missing model1.pdb".format(protein_id))
            continue
        if not _nonempty_file(seq):
            skipped.append("{}: missing seq.fasta".format(protein_id))
            continue

        workdir = output_dir / protein_id
        consensus = workdir / "consensus.tsv"
        if _nonempty_file(consensus) and not force:
            skipped.append("{}: consensus.tsv already exists".format(protein_id))
            continue
        workdir.mkdir(parents=True, exist_ok=True)
        _link_or_copy(model, workdir / "input.pdb", link_inputs)
        _link_or_copy(seq, workdir / "seq.fasta", link_inputs)
        script = _write_starfunc_sbatch(config, protein_id, workdir)
        prepared.append(script)
        if submit:
            submit_script(config, script, workdir)

    return FunctionPreparationResult(output_dir=output_dir, prepared=prepared, skipped=skipped)


def _discover_structure_ids(structure_dir: Path) -> List[str]:
    """Return structure IDs when no decision-tree output is available."""
    return sorted(
        path.name
        for path in structure_dir.iterdir()
        if path.is_dir() and (path / "model1.pdb").exists() and (path / "seq.fasta").exists()
    )


def _write_starfunc_sbatch(config: Dict[str, Any], protein_id: str, workdir: Path) -> Path:
    function_cfg = config.get("function_prediction", {})
    starfunc_sif = _required_config(function_cfg, "starfunc_sif")
    starfunc_database = _required_config(function_cfg, "starfunc_database")
    singularity_module = function_cfg.get("singularity_module", "singularity")
    singularity_command = function_cfg.get("singularity_command", "singularity")
    resources = function_cfg.get("slurm", {})
    body = """#!/bin/bash
set -euo pipefail
{slurm_header}

cd {workdir}
{module_load}

{singularity} run {starfunc_sif} ./ {starfunc_database}
""".format(
        protein_id=protein_id,
        ntasks_per_node=resources.get("ntasks_per_node", 8),
        time=resources.get("time", "10:00:00"),
        mem=resources.get("mem", "100G"),
        slurm_header=slurm_header(
            config,
            job_name="{}_StarFunc".format(protein_id),
            stdout="{}_StarFunc.out".format(protein_id),
            stderr="{}_StarFunc.err".format(protein_id),
            ntasks_per_node=resources.get("ntasks_per_node", 8),
            time=resources.get("time", "10:00:00"),
            mem=resources.get("mem", "100G"),
            partition=resources.get("partition"),
        ),
        workdir=shlex.quote(str(workdir)),
        module_load=module_load_line(config, singularity_module),
        singularity=shlex.quote(singularity_command),
        starfunc_sif=shlex.quote(starfunc_sif),
        starfunc_database=shlex.quote(starfunc_database),
    )
    path = workdir / ("run_starfunc" + script_suffix(config))
    path.write_text(body, encoding="utf-8")
    os.chmod(str(path), 0o755)
    return path


def _nonempty_file(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 1


def _link_or_copy(source: Path, target: Path, use_symlink: bool) -> None:
    if target.exists() or target.is_symlink():
        return
    if use_symlink:
        target.symlink_to(source)
        return
    import shutil

    shutil.copy2(str(source), str(target))


def _required_config(config: Dict[str, Any], key: str) -> str:
    if key not in config or not config[key]:
        raise KeyError("Missing required function_prediction.{} config".format(key))
    return config[key]
