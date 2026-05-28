"""Prepare structure-prediction jobs from decision-tree output."""

from pathlib import Path
from typing import Any, Dict, List, Optional
import os
import shlex

from .config import get_required
from .decisions import Decision, afdb_model_url, read_decision_map
from .execution import module_load_line, script_suffix, slurm_header, submit_script
from .fasta import protein_id_from_header, read_fasta_records, write_single_record


class StructurePreparationResult(object):
    def __init__(self, output_dir: Path, prepared: List[Path], skipped: List[str]) -> None:
        self.output_dir = output_dir
        self.prepared = prepared
        self.skipped = skipped


def prepare_structure_predictions(
    config: Dict[str, Any], submit: bool = False, force: bool = False
) -> StructurePreparationResult:
    """Prepare per-protein structure inputs and Slurm scripts."""
    run_dir = Path(get_required(config, "run.work_dir")).expanduser().resolve()
    input_fasta = Path(get_required(config, "run.input_fasta")).expanduser().resolve()
    structure_cfg = config.get("structure_prediction", {})
    output_dir = Path(
        structure_cfg.get("output_dir", str(run_dir / "structures"))
    ).expanduser().resolve()
    if not input_fasta.is_file():
        raise FileNotFoundError("Missing input FASTA: {}".format(input_fasta))

    output_dir.mkdir(parents=True, exist_ok=True)
    sequences = {
        protein_id_from_header(header): sequence
        for header, sequence in read_fasta_records(input_fasta)
    }
    user_models = structure_cfg.get("user_models", {})
    if isinstance(user_models, dict) and user_models and user_models.get("enabled", True):
        return _prepare_user_model_structures(
            output_dir=output_dir,
            sequences=sequences,
            user_models=user_models,
            force=force,
        )

    decision_file = Path(
        structure_cfg.get(
            "decision_file",
            str(run_dir / "decision_tree_intermediates" / "decide" / "decisions.txt"),
        )
    ).expanduser().resolve()
    enabled_steps = set(structure_cfg.get("enabled_steps", ["GETSEQS", "AFDB"]))

    if not decision_file.is_file():
        raise FileNotFoundError("Missing decision file: {}".format(decision_file))

    decisions = read_decision_map(decision_file)

    prepared = []  # type: List[Path]
    skipped = []  # type: List[str]
    for protein_id, decision in sorted(decisions.items()):
        protein_dir = output_dir / protein_id
        protein_dir.mkdir(parents=True, exist_ok=True)

        sequence = sequences.get(protein_id)
        if sequence is None:
            skipped.append("{}: missing sequence in FASTA".format(protein_id))
            continue
        if "GETSEQS" in enabled_steps or not (protein_dir / "seq.fasta").exists():
            write_single_record(protein_dir / "seq.fasta", protein_id, sequence)

        model_path = protein_dir / "model1.pdb"
        if model_path.exists() and model_path.stat().st_size > 1 and not force:
            skipped.append("{}: model1.pdb already exists".format(protein_id))
            continue

        script_path = _prepare_one_structure_job(
            config=config,
            decision=decision,
            protein_dir=protein_dir,
            enabled_steps=enabled_steps,
        )
        if script_path is None:
            skipped.append(
                "{}: no enabled structure action for decision {}".format(
                    protein_id, decision.method
                )
            )
            continue
        prepared.append(script_path)
        if submit and script_path.suffix in (".sbatch", ".sh"):
            submit_script(config, script_path, protein_dir)
    return StructurePreparationResult(output_dir=output_dir, prepared=prepared, skipped=skipped)


def _prepare_user_model_structures(
    output_dir: Path,
    sequences: Dict[str, str],
    user_models: Dict[str, Any],
    force: bool,
) -> StructurePreparationResult:
    """Stage user-provided structures into the standard StarFunc layout."""
    if not sequences:
        raise ValueError("Input FASTA did not contain any records")
    model_map = _read_user_model_mapping(user_models)
    model_dir = None  # type: Optional[Path]
    if user_models.get("model_dir"):
        model_dir = Path(user_models["model_dir"]).expanduser().resolve()
    pattern = user_models.get("pattern", "{protein_id}.pdb")
    link_models = bool(user_models.get("link_models", user_models.get("link_inputs", True)))

    prepared = []  # type: List[Path]
    skipped = []  # type: List[str]
    for protein_id, sequence in sorted(sequences.items()):
        protein_dir = output_dir / protein_id
        protein_dir.mkdir(parents=True, exist_ok=True)
        write_single_record(protein_dir / "seq.fasta", protein_id, sequence)

        source = _resolve_user_model_source(protein_id, model_map, model_dir, pattern)
        if source is None:
            skipped.append("{}: missing user model mapping".format(protein_id))
            continue
        if not source.is_file() or source.stat().st_size < 1:
            skipped.append("{}: missing user model {}".format(protein_id, source))
            continue

        target = protein_dir / "model1.pdb"
        if target.exists() or target.is_symlink():
            if not force:
                skipped.append("{}: model1.pdb already exists".format(protein_id))
                continue
            target.unlink()
        _link_or_copy(source, target, link_models)
        marker = protein_dir / "use_user_model.done"
        marker.write_text(str(source) + "\n", encoding="utf-8")
        prepared.append(marker)

    return StructurePreparationResult(output_dir=output_dir, prepared=prepared, skipped=skipped)


def _read_user_model_mapping(user_models: Dict[str, Any]) -> Dict[str, Path]:
    mapping_file = user_models.get("mapping_file")
    if not mapping_file:
        return {}
    path = Path(mapping_file).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError("Missing user model mapping file: {}".format(path))
    mapping = {}  # type: Dict[str, Path]
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(
                    "Expected at least protein_id and model_path in {} line {}".format(
                        path, line_number
                    )
                )
            if parts[0].lower() in ("protein_id", "id"):
                continue
            model_path = Path(parts[1]).expanduser()
            if not model_path.is_absolute():
                model_path = path.parent / model_path
            mapping[parts[0]] = model_path.resolve()
    return mapping


def _resolve_user_model_source(
    protein_id: str,
    model_map: Dict[str, Path],
    model_dir: Optional[Path],
    pattern: str,
) -> Optional[Path]:
    if protein_id in model_map:
        return model_map[protein_id]
    if model_dir is None:
        return None
    relative = Path(pattern.format(protein_id=protein_id))
    if relative.is_absolute():
        return relative.expanduser().resolve()
    return (model_dir / relative).expanduser().resolve()


def _prepare_one_structure_job(
    config: Dict[str, Any],
    decision: Decision,
    protein_dir: Path,
    enabled_steps: set,
) -> Optional[Path]:
    if decision.method == "AFDB" and "AFDB" in enabled_steps:
        return _prepare_afdb_model(config, decision, protein_dir)
    if decision.method == "MODELLER" and "MODELLER" in enabled_steps:
        template_url = _afdb_url_for_decision(decision)
        return _write_modeler_sbatch(config, decision.protein_id, protein_dir, template_url)
    if decision.method == "LOMETS" and "LOMETS" in enabled_steps:
        return _write_lomets_sbatch(config, decision.protein_id, protein_dir)
    if decision.method == "DITASSER" and "DITASSER" in enabled_steps:
        return _write_ditasser_sbatch(config, decision.protein_id, protein_dir)
    if "DMFOLD" in enabled_steps:
        return _write_dmfold_sbatch(config, decision.protein_id, protein_dir)
    return None


def _prepare_afdb_model(
    config: Dict[str, Any],
    decision: Decision,
    protein_dir: Path,
) -> Path:
    url = _afdb_url_for_decision(decision)
    script = protein_dir / ("fetch_afdb" + script_suffix(config))
    script.write_text(
        _afdb_fetch_script(
            config,
            decision.protein_id,
            protein_dir,
            url,
            "model1.pdb",
        ),
        encoding="utf-8",
    )
    return script


def _link_or_copy(source: Path, target: Path, use_symlink: bool) -> None:
    if target.exists() or target.is_symlink():
        return
    if use_symlink:
        target.symlink_to(source)
        return
    import shutil

    shutil.copy2(str(source), str(target))


def _write_ditasser_sbatch(config: Dict[str, Any], protein_id: str, protein_dir: Path) -> Path:
    tools = _structure_tools(config)
    body = """#!/bin/bash
set -euo pipefail
{slurm_header}

cd {workdir}
{module_load}

{singularity} exec -B {scratch_bind}:/tmp -B {pkgdir}:/D-ITASSER-2.0:ro -B {itlib}:/ITLIB:ro -B {workdir}:/datadir {sif} perl /I-TASSERmod/runI-TASSER.pl -libdir /ITLIB -seqname {protein_id} -datadir /datadir -pkgdir /D-ITASSER-2.0 -outdir /datadir -runstyle gnuparallel -nmodel 1 -itmode "DIT-AF2" -msapipe "DeepMSA2" > run.log 2> run.err
""".format(
        protein_id=protein_id,
        slurm_header=slurm_header(
            config,
            job_name="{}_itas".format(protein_id),
            stdout="{}_itasser.out".format(protein_id),
            stderr="{}_itasser.err".format(protein_id),
            ntasks_per_node=8,
            time="72:00:00",
            mem="100G",
        ),
        workdir=shlex.quote(str(protein_dir)),
        scratch_bind=shlex.quote(tools["scratch_bind"]),
        pkgdir=shlex.quote(tools["ditasser_pkgdir"]),
        itlib=shlex.quote(tools["itlib_dir"]),
        sif=shlex.quote(tools["seq2fun_sif"]),
        singularity=shlex.quote(tools["singularity_command"]),
        module_load=module_load_line(config, tools["singularity_module"]),
    )
    return _write_script(protein_dir / ("run_ditasser" + script_suffix(config)), body)


def _write_lomets_sbatch(config: Dict[str, Any], protein_id: str, protein_dir: Path) -> Path:
    tools = _structure_tools(config)
    body = """#!/bin/bash
set -euo pipefail
{slurm_header}

cd {workdir}
{module_load}

{singularity} exec -B {scratch_bind}:/tmp -B {pkgdir}:/D-ITASSER-2.0:ro -B {itlib}:/ITLIB:ro -B {workdir}:/datadir {sif} perl /I-TASSERmod/runLOMETS.pl -libdir /ITLIB -seqname {protein_id} -datadir /datadir -pkgdir /D-ITASSER-2.0 -outdir /datadir -runstyle gnuparallel
{singularity} exec -B {scratch_bind}:/tmp -B {workdir}:{workdir} -B {pkgdir}:/D-ITASSER-2.0:ro {sif} perl /I-TASSERmod/modeler/initdat2model_noss.pl {workdir}/seq.fasta {workdir}/initall.dat
""".format(
        protein_id=protein_id,
        slurm_header=slurm_header(
            config,
            job_name="{}_lomets".format(protein_id),
            stdout="{}_lomets.out".format(protein_id),
            stderr="{}_lomets.err".format(protein_id),
            ntasks_per_node=8,
            time="72:00:00",
            mem="50G",
        ),
        workdir=shlex.quote(str(protein_dir)),
        scratch_bind=shlex.quote(tools["scratch_bind"]),
        pkgdir=shlex.quote(tools["ditasser_pkgdir"]),
        itlib=shlex.quote(tools["itlib_dir"]),
        sif=shlex.quote(tools["seq2fun_sif"]),
        singularity=shlex.quote(tools["singularity_command"]),
        module_load=module_load_line(config, tools["singularity_module"]),
    )
    return _write_script(protein_dir / ("run_lomets" + script_suffix(config)), body)


def _write_modeler_sbatch(
    config: Dict[str, Any],
    protein_id: str,
    protein_dir: Path,
    template_url: Optional[str] = None,
) -> Path:
    tools = _structure_tools(config)
    fetch_template = ""
    timeout = _model_download_timeout(config)
    if template_url:
        fetch_template = """
if [ ! -s template.pdb ]; then
    python3 - <<'PY'
from pathlib import Path
import shutil
import urllib.parse
import urllib.request

url = {template_url!r}
timeout = {timeout!r}
suffix = Path(urllib.parse.urlparse(url).path).suffix.lower()
download = "template" + (suffix if suffix else ".pdb")
with urllib.request.urlopen(url, timeout=timeout) as response, open(download, "wb") as handle:
    shutil.copyfileobj(response, handle)
if suffix in (".cif", ".bcif"):
    import gemmi
    structure = gemmi.read_structure(download)
    structure.write_pdb("template.pdb")
else:
    Path(download).replace("template.pdb")
PY
fi
""".format(template_url=template_url, timeout=timeout)
    body = """#!/bin/bash
set -euo pipefail
{slurm_header}

cd {workdir}
{module_load}
{fetch_template}

{singularity} exec -B {scratch_bind}:/tmp -B {workdir}:{workdir} -B {pkgdir}:/D-ITASSER-2.0:ro {sif} perl /I-TASSERmod/modeler/modeller_sw.pl {workdir}/seq.fasta {workdir}/template.pdb {workdir}/model1.pdb 0
""".format(
        protein_id=protein_id,
        slurm_header=slurm_header(
            config,
            job_name="{}_modeller".format(protein_id),
            stdout="{}_modeller.out".format(protein_id),
            stderr="{}_modeller.err".format(protein_id),
            ntasks=1,
            time="16:00:00",
            mem="10G",
        ),
        workdir=shlex.quote(str(protein_dir)),
        scratch_bind=shlex.quote(tools["scratch_bind"]),
        pkgdir=shlex.quote(tools["ditasser_pkgdir"]),
        sif=shlex.quote(tools["seq2fun_sif"]),
        singularity=shlex.quote(tools["singularity_command"]),
        module_load=module_load_line(config, tools["singularity_module"]),
        fetch_template=fetch_template,
    )
    return _write_script(protein_dir / ("run_modeler" + script_suffix(config)), body)


def _write_dmfold_sbatch(config: Dict[str, Any], protein_id: str, protein_dir: Path) -> Path:
    tools = _structure_tools(config)
    conda = shlex.quote(tools["conda_executable"])
    env = shlex.quote(tools["dmfold_conda_env"])
    base = shlex.quote(tools["dmfold_base"])
    workdir = shlex.quote(str(protein_dir))
    body = """#!/bin/bash
set -euo pipefail
{slurm_header}

cd {workdir}

{conda} run -p {env} python3 {base}/Parse_sequence.py -o={workdir}
{conda} run -p {env} python3 {base}/Run_DeepMSA2.py -o={workdir}
{conda} run -p {env} python3 {base}/MSA_combination.py -o={workdir}
{conda} run -p {env} python3 {base}/Run_DMFold.py -o={workdir}
""".format(
        protein_id=protein_id,
        slurm_header=slurm_header(
            config,
            job_name="{}_dmfold".format(protein_id),
            stdout="{}_dmfold.out".format(protein_id),
            stderr="{}_dmfold.err".format(protein_id),
            ntasks_per_node=8,
            time="72:00:00",
            mem="50G",
        ),
        workdir=workdir,
        conda=conda,
        env=env,
        base=base,
    )
    return _write_script(protein_dir / ("run_dmfold" + script_suffix(config)), body)


def _afdb_fetch_script(
    config: Dict[str, Any],
    protein_id: str,
    protein_dir: Path,
    model_url: str,
    output_name: str,
) -> str:
    timeout = _model_download_timeout(config)
    return """#!/bin/bash
set -euo pipefail
{slurm_header}

cd {workdir}

python3 - <<'PY'
from pathlib import Path
import shutil
import urllib.parse
import urllib.request

url = {model_url!r}
timeout = {timeout!r}
output = Path({output_name!r})
suffix = Path(urllib.parse.urlparse(url).path).suffix.lower()
download = output.with_suffix(suffix if suffix else ".pdb")
with urllib.request.urlopen(url, timeout=timeout) as response, open(download, "wb") as handle:
    shutil.copyfileobj(response, handle)
if suffix in (".cif", ".bcif"):
    import gemmi
    structure = gemmi.read_structure(str(download))
    structure.write_pdb(str(output))
else:
    download.replace(output)
PY
""".format(
        protein_id=protein_id,
        slurm_header=slurm_header(
            config,
            job_name="{}_afdb".format(protein_id),
            stdout="{}_afdb.out".format(protein_id),
            stderr="{}_afdb.err".format(protein_id),
            ntasks=1,
            time="01:00:00",
            mem="1G",
        ),
        workdir=shlex.quote(str(protein_dir)),
        model_url=model_url,
        timeout=timeout,
        output_name=output_name,
    )


def _afdb_url_for_decision(decision: Decision) -> str:
    url = afdb_model_url(decision.match)
    if url:
        return url
    model_id = _afdb_model_id(decision.match)
    if model_id is None:
        raise ValueError("Decision does not contain an AFDB match: {}".format(decision))
    return "https://alphafold.ebi.ac.uk/files/{}-model_v6.cif".format(model_id)


def _afdb_model_id(match: str) -> Optional[str]:
    if not match or not match.startswith("AFDB:"):
        return None
    return _strip_afdb_model_suffix(match.split(":", 1)[1].split("|", 1)[0])


def _strip_afdb_model_suffix(match: str) -> str:
    for suffix in (".pdb", ".cif", ".bcif"):
        if match.endswith(suffix):
            match = match[: -len(suffix)]
            break
    marker = "-model_v"
    if marker in match:
        return match.split(marker, 1)[0]
    return match


def _model_download_timeout(config: Dict[str, Any]) -> float:
    structure_cfg = config.get("structure_prediction", {})
    decision_cfg = config.get("decision_tree", {})
    return float(structure_cfg.get("model_download_timeout", decision_cfg.get("afdb_api_timeout", 30)))


def _write_script(path: Path, body: str) -> Path:
    path.write_text(body, encoding="utf-8")
    os.chmod(str(path), 0o755)
    return path


def _structure_tools(config: Dict[str, Any]) -> Dict[str, str]:
    cfg = config.get("structure_prediction", {}).get("tools", {})
    return {
        "seq2fun_sif": _required_tool(cfg, "seq2fun_sif"),
        "ditasser_pkgdir": _required_tool(cfg, "ditasser_pkgdir"),
        "itlib_dir": _required_tool(cfg, "itlib_dir"),
        "scratch_bind": _required_tool(cfg, "scratch_bind"),
        "dmfold_base": cfg.get("dmfold_base", ""),
        "dmfold_conda_env": cfg.get("dmfold_conda_env", ""),
        "conda_executable": config.get("decision_tree", {}).get("conda_executable", ""),
        "singularity_module": cfg.get("singularity_module", "singularity"),
        "singularity_command": cfg.get("singularity_command", "singularity"),
    }


def _required_tool(cfg: Dict[str, Any], key: str) -> str:
    if key not in cfg or not cfg[key]:
        raise KeyError("Missing required structure_prediction.tools.{} config".format(key))
    return cfg[key]
