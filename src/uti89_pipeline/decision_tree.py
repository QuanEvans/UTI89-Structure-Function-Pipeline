"""Prepare and run the decision-tree stage."""

import json
import os
import shlex
import subprocess
import filecmp
import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional

from .config import get_required, write_yaml
from .execution import snakemake_profile


class DecisionTreeRun:
    def __init__(
        self,
        run_dir: Path,
        work_dir: Path,
        pipeline_files: Path,
        input_fasta: Path,
        decisions: Path,
        command: List[str],
    ) -> None:
        self.run_dir = run_dir
        self.work_dir = work_dir
        self.pipeline_files = pipeline_files
        self.input_fasta = input_fasta
        self.decisions = decisions
        self.command = command


def prepare_decision_tree_run(config: Dict[str, Any]) -> DecisionTreeRun:
    """Create an isolated Snakemake work directory for one decision-tree run."""
    run_dir = Path(get_required(config, "run.work_dir")).expanduser().resolve()
    input_fasta = Path(get_required(config, "run.input_fasta")).expanduser().resolve()
    source_dir = Path(get_required(config, "decision_tree.source_dir")).resolve()
    conda_env = Path(get_required(config, "decision_tree.conda_env")).expanduser().resolve()

    _require_file(input_fasta, "input FASTA")
    _require_dir(source_dir, "decision-tree source directory")

    work_dir = run_dir / "decision_tree"
    pipeline_files = run_dir / "decision_tree_intermediates"
    work_dir.mkdir(parents=True, exist_ok=True)
    pipeline_files.mkdir(parents=True, exist_ok=True)
    staged_fasta = pipeline_files / "input.fasta"
    _stage_input_fasta(input_fasta, staged_fasta)

    _link_source_entry(source_dir / "snakefile", work_dir / "snakefile")
    _link_source_entry(source_dir / "workflow", work_dir / "workflow")
    _link_source_entry(source_dir / "scripts", work_dir / "scripts")
    if (source_dir / "environment.yaml").exists():
        _link_source_entry(source_dir / "environment.yaml", work_dir / "environment.yaml")

    profile_dir = work_dir / ".smk_profile"
    profile_dir.mkdir(exist_ok=True)
    cluster_status = source_dir / ".smk_profile" / "cluster_status.py"
    if cluster_status.exists():
        _link_source_entry(cluster_status, profile_dir / "cluster_status.py")

    decision_config = _build_decision_tree_config(config, work_dir, pipeline_files, conda_env)
    write_yaml(work_dir / "config.yaml", decision_config)
    _write_snakemake_profile(config, profile_dir, cluster_status.exists())

    python_bin = conda_env / "bin" / "python"
    command = [
        str(python_bin),
        "-m",
        "snakemake",
        "--snakefile",
        str(work_dir / "workflow" / "1_decide.smk"),
        "--profile",
        ".smk_profile",
        "decideall.done",
    ]

    decisions = pipeline_files / "decide" / "decisions.txt"
    run = DecisionTreeRun(
        run_dir=run_dir,
        work_dir=work_dir,
        pipeline_files=pipeline_files,
        input_fasta=input_fasta,
        decisions=decisions,
        command=command,
    )
    _write_metadata(run, source_dir)
    return run


def run_decision_tree(
    run: DecisionTreeRun,
    dry_run: bool = False,
    unlock: bool = False,
    extra_args: Optional[List[str]] = None,
) -> int:
    """Run Snakemake for a prepared decision-tree run."""
    command = list(run.command)
    if dry_run:
        command.append("-n")
    if unlock:
        command.append("--unlock")
    if extra_args:
        command.extend(extra_args)

    python_bin = Path(command[0])
    _require_file(python_bin, "decision-tree Python executable")
    env = os.environ.copy()
    cache_dir = run.run_dir / ".cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    env.setdefault("XDG_CACHE_HOME", str(cache_dir))
    completed = subprocess.run(command, cwd=run.work_dir, env=env, check=False)
    return completed.returncode


def format_command(command: List[str], cwd: Path) -> str:
    """Return a shell-ready command display string."""
    return "cd {} && {}".format(
        shlex.quote(str(cwd)), " ".join(shlex.quote(part) for part in command)
    )


def _build_decision_tree_config(
    config: Dict[str, Any],
    work_dir: Path,
    pipeline_files: Path,
    conda_env: Path,
) -> Dict[str, Any]:
    decision_cfg = config["decision_tree"]
    return {
        "PDB_MIN_SEQ_ID": decision_cfg.get("pdb_min_seq_id", 0.95),
        "AF_MIN_SEQ_ID": decision_cfg.get("af_min_seq_id", 0.9),
        "AF_MIN_PLDDT": decision_cfg.get("af_min_plddt", 70),
        "EXACT_SEQ_ID_DEF": decision_cfg.get("exact_seq_id_def", 1.0),
        "PATH_INFASTA": str(pipeline_files / "input.fasta"),
        "BASE_PATH": str(work_dir),
        "PATH_SUBWORKFLOW_0": str(work_dir / "workflow" / "0_download_fasta_build_dmnd.smk"),
        "PATH_SUBWORKFLOW_1": str(work_dir / "workflow" / "1_decide.smk"),
        "PATH_SUBWORKFLOW_2": str(work_dir / "workflow" / "2_predict.smk"),
        "PATH_SCRIPTS": str(work_dir / "scripts"),
        "PATH_PIPELINE_FILES": str(pipeline_files),
        "PATH_DMND_LIB": str(Path(get_required(config, "decision_tree.dmnd_lib")).expanduser()),
        "AFDB_FOLDCOMP_PATH": str(
            Path(get_required(config, "decision_tree.afdb_foldcomp")).expanduser()
        ),
        "PATH_FOLDCOMP_SUBDB": str(pipeline_files / "foldcomp" / "afdb_matches_subdb"),
        "PATH_FASTA_LIB": str(Path(get_required(config, "decision_tree.fasta_lib")).expanduser()),
        "PATH_CONDA": str(conda_env),
        "URL_AFDB_FASTA": decision_cfg.get(
            "url_afdb_fasta", "https://ftp.ebi.ac.uk/pub/databases/alphafold/sequences.fasta"
        ),
        "URL_PDB_FASTA": decision_cfg.get(
            "url_pdb_fasta", "https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt"
        ),
        "URL_PDB_SEQ_INFO": decision_cfg.get(
            "url_pdb_seq_info",
            "https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_entry_type.txt",
        ),
    }


def _write_snakemake_profile(
    config: Dict[str, Any], profile_dir: Path, include_cluster_status: bool
) -> None:
    write_yaml(profile_dir / "config.yaml", snakemake_profile(config, include_cluster_status))


def _write_metadata(run: DecisionTreeRun, source_dir: Path) -> None:
    metadata = {
        "stage": "decision_tree",
        "source_dir": str(source_dir),
        "run_dir": str(run.run_dir),
        "work_dir": str(run.work_dir),
        "pipeline_files": str(run.pipeline_files),
        "input_fasta": str(run.input_fasta),
        "decisions": str(run.decisions),
        "command": run.command,
        "command_display": format_command(run.command, run.work_dir),
    }
    (run.run_dir / "decision_tree_run.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )


def _link_source_entry(source: Path, target: Path) -> None:
    if not source.exists():
        raise FileNotFoundError(f"Required source entry does not exist: {source}")
    if target.exists() or target.is_symlink():
        if target.is_symlink() and Path(os.readlink(target)) == source:
            return
        raise FileExistsError(f"Refusing to overwrite existing path: {target}")
    target.symlink_to(source, target_is_directory=source.is_dir())


def _stage_input_fasta(input_fasta: Path, staged_fasta: Path) -> None:
    """Copy input FASTA once and refuse silent input changes for a run dir."""
    if not staged_fasta.exists():
        shutil.copy2(str(input_fasta), str(staged_fasta))
        return
    if not filecmp.cmp(str(input_fasta), str(staged_fasta), shallow=False):
        raise ValueError(
            "Input FASTA differs from the staged FASTA for this run directory. "
            "Use a fresh run.work_dir or pass a different --work-dir."
        )


def _require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {label}: {path}")


def _require_dir(path: Path, label: str) -> None:
    if not path.is_dir():
        raise NotADirectoryError(f"Missing {label}: {path}")
