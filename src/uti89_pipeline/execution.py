"""Execution backend helpers."""

import subprocess
from pathlib import Path
from typing import Any, Dict, Optional


def backend(config: Dict[str, Any]) -> str:
    """Return the configured execution backend."""
    value = config.get("execution", {}).get("backend")
    if value:
        return str(value).lower()
    return "slurm" if "slurm" in config else "local"


def is_slurm(config: Dict[str, Any]) -> bool:
    return backend(config) == "slurm"


def local_cores(config: Dict[str, Any]) -> int:
    return int(config.get("execution", {}).get("local_cores", 1))


def module_load_line(config: Dict[str, Any], module_name: Optional[str]) -> str:
    """Return a module-load line only when requested for the backend."""
    if not module_name:
        return ""
    execution_cfg = config.get("execution", {})
    use_modules = execution_cfg.get("use_modules")
    if use_modules is None:
        use_modules = is_slurm(config)
    if not use_modules:
        return ""
    return "module load {}\n".format(module_name)


def submit_script(config: Dict[str, Any], script: Path, cwd: Path) -> None:
    """Submit or run a generated script according to the configured backend."""
    if is_slurm(config):
        command = ["sbatch", str(script)]
    else:
        command = ["bash", str(script)]
    subprocess.check_call(command, cwd=str(cwd))


def slurm_header(
    config: Dict[str, Any],
    job_name: str,
    stdout: str,
    stderr: str,
    nodes: int = 1,
    ntasks: Optional[int] = None,
    ntasks_per_node: Optional[int] = None,
    time: str = "01:00:00",
    mem: str = "1G",
    partition: Optional[str] = None,
) -> str:
    """Return SBATCH directives, or an empty string for local execution."""
    if not is_slurm(config):
        return ""
    slurm = _slurm_config(config)
    lines = [
        "#SBATCH --job-name={}".format(job_name),
        "#SBATCH --account={}".format(slurm["account"]),
        "#SBATCH --partition={}".format(partition or slurm["partition"]),
        "#SBATCH --output={}".format(stdout),
        "#SBATCH --error={}".format(stderr),
        "#SBATCH --nodes={}".format(nodes),
    ]
    if ntasks is not None:
        lines.append("#SBATCH --ntasks={}".format(ntasks))
    if ntasks_per_node is not None:
        lines.append("#SBATCH --ntasks-per-node={}".format(ntasks_per_node))
    lines.extend(
        [
            "#SBATCH --time={}".format(time),
            "#SBATCH --mem={}".format(mem),
        ]
    )
    return "\n".join(lines) + "\n"


def script_suffix(config: Dict[str, Any]) -> str:
    return ".sbatch" if is_slurm(config) else ".sh"


def _slurm_config(config: Dict[str, Any]) -> Dict[str, Any]:
    slurm = config.get("slurm")
    if not isinstance(slurm, dict) or not slurm.get("account") or not slurm.get("partition"):
        raise KeyError(
            "execution.backend is slurm, so config must include slurm.account and slurm.partition"
        )
    return slurm
