"""Prepare and run the redesigned native decision-tree stage."""

import csv
import filecmp
import json
import shlex
import shutil
import subprocess
import sys
import urllib.request
from urllib.parse import quote, urlencode
from pathlib import Path
from typing import Any, Dict, Iterable, List, NamedTuple, Optional, Tuple

from .config import get_required
from .fasta import protein_id_from_header, read_fasta_records


class SequenceRecord(NamedTuple):
    protein_id: str
    header: str
    sequence: str


class Match(NamedTuple):
    protein_id: str
    match: str
    evalue: Optional[float]
    seq_id: float
    query_coverage: float
    target_coverage: float
    query_length: int
    target_length: int
    align_length: int
    plddt: Optional[float] = None
    pdb_path: Optional[str] = None


class DecisionTreeRun:
    def __init__(
        self,
        config: Dict[str, Any],
        run_dir: Path,
        work_dir: Path,
        pipeline_files: Path,
        input_fasta: Path,
        decisions: Path,
        command: List[str],
    ) -> None:
        self.config = config
        self.run_dir = run_dir
        self.work_dir = work_dir
        self.pipeline_files = pipeline_files
        self.input_fasta = input_fasta
        self.decisions = decisions
        self.command = command


def prepare_decision_tree_run(config: Dict[str, Any]) -> DecisionTreeRun:
    """Create the native decision-tree run directory and metadata."""
    run_dir = Path(get_required(config, "run.work_dir")).expanduser().resolve()
    input_fasta = Path(get_required(config, "run.input_fasta")).expanduser().resolve()
    _require_file(input_fasta, "input FASTA")

    work_dir = run_dir / "decision_tree"
    pipeline_files = run_dir / "decision_tree_intermediates"
    decide_dir = pipeline_files / "decide"
    reports_dir = pipeline_files / "reports"
    work_dir.mkdir(parents=True, exist_ok=True)
    decide_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    staged_fasta = pipeline_files / "input.fasta"
    _stage_input_fasta(input_fasta, staged_fasta)

    decisions = decide_dir / "decisions.txt"
    command = [
        sys.executable,
        "scripts/run_decision_tree.py",
        "--config",
        "<config>",
    ]
    run = DecisionTreeRun(
        config=config,
        run_dir=run_dir,
        work_dir=work_dir,
        pipeline_files=pipeline_files,
        input_fasta=input_fasta,
        decisions=decisions,
        command=command,
    )
    _write_metadata(run)
    return run


def run_decision_tree(
    run: DecisionTreeRun,
    dry_run: bool = False,
    unlock: bool = False,
    extra_args: Optional[List[str]] = None,
) -> int:
    """Run the native decision tree and write decisions.txt."""
    if unlock:
        return 0
    if extra_args:
        raise ValueError("Native decision tree does not accept extra command-line arguments")
    if dry_run:
        _write_plan(run)
        return 0
    _write_decisions(run)
    return 0


def format_command(command: List[str], cwd: Path) -> str:
    """Return a shell-ready command display string."""
    return "cd {} && {}".format(
        shlex.quote(str(cwd)), " ".join(shlex.quote(part) for part in command)
    )


def choose_structure_method(
    sequence: str,
    af_match: Optional[Match],
    pdb_match: Optional[Match],
    config: Dict[str, Any],
) -> Tuple[str, str, str, str]:
    """Return method, match, score, reason for one protein."""
    decision_cfg = config.get("decision_tree", {})
    structure_cfg = config.get("structure_prediction", {})
    enabled_steps = set(structure_cfg.get("enabled_steps", ["GETSEQS", "AFDB", "MODELLER", "LOMETS", "DITASSER"]))

    af_exact_seq_id = float(decision_cfg.get("af_exact_seq_id", 1.0))
    af_min_seq_id = float(decision_cfg.get("af_min_seq_id", 0.9))
    af_min_query_coverage = float(decision_cfg.get("af_min_query_coverage", 0.8))
    af_min_target_coverage = float(decision_cfg.get("af_min_target_coverage", 0.8))
    af_min_plddt = float(decision_cfg.get("af_min_plddt", 70.0))
    pdb_min_seq_id = float(decision_cfg.get("pdb_min_seq_id", 0.95))
    pdb_min_query_coverage = float(decision_cfg.get("pdb_min_query_coverage", 0.8))
    long_sequence_length = int(decision_cfg.get("long_sequence_length", 900))
    long_sequence_method = str(decision_cfg.get("long_sequence_method", "DMFOLD")).upper()

    if af_match is not None and _has_good_plddt(af_match, af_min_plddt):
        if (
            af_match.seq_id >= af_exact_seq_id
            and af_match.query_coverage >= 0.99
            and af_match.target_coverage >= 0.99
            and "AFDB" in enabled_steps
        ):
            return "AFDB", af_match.match, _score(af_match.seq_id), "afdb_exact_high_confidence"

    if (
        pdb_match is not None
        and pdb_match.seq_id >= pdb_min_seq_id
        and pdb_match.query_coverage >= pdb_min_query_coverage
        and "LOMETS" in enabled_steps
    ):
        return "LOMETS", pdb_match.match, _score(pdb_match.seq_id), "pdb_high_identity"

    if af_match is not None and _has_good_plddt(af_match, af_min_plddt):
        if (
            af_match.seq_id >= af_min_seq_id
            and af_match.query_coverage >= af_min_query_coverage
            and af_match.target_coverage >= af_min_target_coverage
            and "MODELLER" in enabled_steps
        ):
            return "MODELLER", af_match.match, _score(af_match.seq_id), "afdb_template_high_confidence"

    if len(sequence) >= long_sequence_length and long_sequence_method in enabled_steps:
        return long_sequence_method, "", "0", "long_sequence_fallback"
    if "DITASSER" in enabled_steps:
        return "DITASSER", "", "0", "no_high_confidence_template"
    if "DMFOLD" in enabled_steps:
        return "DMFOLD", "", "0", "fallback"
    return "NEEDS_USER_MODEL", "", "0", "no_enabled_structure_method"


def _write_decisions(run: DecisionTreeRun) -> None:
    records = _read_sequence_records(run.pipeline_files / "input.fasta")
    pdb_matches = _load_pdb_matches(run, records)
    af_matches = _load_afdb_matches(run, records)

    run.decisions.parent.mkdir(parents=True, exist_ok=True)
    with run.decisions.open("w", encoding="utf-8") as handle:
        for record in records:
            af_match = af_matches.get(record.protein_id)
            pdb_match = pdb_matches.get(record.protein_id)
            method, match, score, reason = choose_structure_method(
                record.sequence,
                af_match,
                pdb_match,
                run.config,
            )
            if method in {"AFDB", "MODELLER"} and af_match is not None:
                match = _decision_match(af_match)
            handle.write(
                "{}\t{}\t{}\t{}\t{}\n".format(
                    record.protein_id, method, match, score, reason
                )
            )


def _write_plan(run: DecisionTreeRun) -> None:
    plan = {
        "stage": "decision_tree",
        "mode": "native",
        "input_fasta": str(run.pipeline_files / "input.fasta"),
        "decisions": str(run.decisions),
        "pdb_search": _pdb_search_plan(run.config),
        "afdb_search": _afdb_search_plan(run.config),
    }
    (run.work_dir / "decision_tree_plan.json").write_text(
        json.dumps(plan, indent=2) + "\n", encoding="utf-8"
    )


def _read_sequence_records(fasta: Path) -> List[SequenceRecord]:
    return [
        SequenceRecord(protein_id_from_header(header), header, sequence)
        for header, sequence in read_fasta_records(fasta)
    ]


def _load_pdb_matches(run: DecisionTreeRun, records: List[SequenceRecord]) -> Dict[str, Match]:
    decision_cfg = run.config.get("decision_tree", {})
    hits_tsv = decision_cfg.get("pdb_hits_tsv")
    if hits_tsv:
        return _read_hits_tsv(Path(hits_tsv).expanduser(), source="PDB")

    pdb_dmnd = _pdb_dmnd_path(decision_cfg)
    if not pdb_dmnd:
        return {}
    output = run.pipeline_files / "reports" / "pdb_hits.tsv"
    if not output.exists():
        _run_diamond(
            query=run.pipeline_files / "input.fasta",
            db=pdb_dmnd,
            output=output,
            threads=int(decision_cfg.get("pdb_threads", decision_cfg.get("threads", 4))),
            diamond_executable=decision_cfg.get("diamond_executable"),
        )
    return _read_hits_tsv(output, source="PDB")


def _load_afdb_matches(run: DecisionTreeRun, records: List[SequenceRecord]) -> Dict[str, Match]:
    decision_cfg = run.config.get("decision_tree", {})
    hits_tsv = decision_cfg.get("afdb_hits_tsv")
    if hits_tsv:
        matches = _read_hits_tsv(Path(hits_tsv).expanduser(), source="AFDB")
    else:
        matches = {}

    use_api = bool(decision_cfg.get("use_afdb_api", True))
    use_sequence_search = bool(decision_cfg.get("afdb_sequence_search", True))
    for record in records:
        if record.protein_id in matches:
            continue
        match = _afdb_api_match(record, decision_cfg) if use_api else None
        if match is None and use_api and use_sequence_search:
            match = _afdb_sequence_match(record, decision_cfg)
        if match is not None:
            matches[record.protein_id] = match
    return matches


def _afdb_api_match(record: SequenceRecord, decision_cfg: Dict[str, Any]) -> Optional[Match]:
    if not _looks_like_uniprot_accession(record.protein_id):
        return None
    base_url = str(decision_cfg.get("afdb_api_base_url", "https://alphafold.ebi.ac.uk/api/prediction")).rstrip("/")
    url = "{}/{}".format(base_url, record.protein_id)
    try:
        with urllib.request.urlopen(url, timeout=float(decision_cfg.get("afdb_api_timeout", 30))) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except Exception:
        return None
    if not payload:
        return None
    entry = payload[0] if isinstance(payload, list) else payload
    target_sequence = entry.get("sequence") or entry.get("uniprotSequence")
    if target_sequence and _clean_sequence(target_sequence) != _clean_sequence(record.sequence):
        return None
    af_id = entry.get("entryId") or "AF-{}-F1".format(record.protein_id)
    plddt = _first_float(entry, ("globalMetricValue", "averagePlddt", "plddt", "confidenceScore"))
    model_url = _first_present(entry, ("cifUrl", "bcifUrl", "pdbUrl"))
    if plddt is None and not model_url:
        return None
    match_name = "AFDB:{}".format(af_id)
    return Match(
        protein_id=record.protein_id,
        match=match_name,
        evalue=0.0,
        seq_id=1.0,
        query_coverage=1.0,
        target_coverage=1.0,
        query_length=len(record.sequence),
        target_length=len(record.sequence),
        align_length=len(record.sequence),
        plddt=plddt,
        pdb_path=model_url,
    )


def _afdb_sequence_match(record: SequenceRecord, decision_cfg: Dict[str, Any]) -> Optional[Match]:
    sequence = _clean_sequence(record.sequence)
    if len(sequence) < int(decision_cfg.get("afdb_sequence_min_length", 20)):
        return None
    if any(residue not in "ACDEFGHIKLMNPQRSTVWY" for residue in sequence):
        return None
    base_url = str(
        decision_cfg.get(
            "afdb_sequence_api_url",
            "https://alphafold.ebi.ac.uk/api/sequence/summary",
        )
    ).rstrip("?")
    params = urlencode(
        {
            "id": sequence,
            "type": "sequence",
            "rows": str(int(decision_cfg.get("afdb_sequence_rows", 20))),
        }
    )
    url = "{}?{}".format(base_url, params)
    try:
        with urllib.request.urlopen(url, timeout=float(decision_cfg.get("afdb_api_timeout", 30))) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except Exception:
        return None
    structures = payload.get("structures") if isinstance(payload, dict) else None
    if not structures:
        return None
    best = None  # type: Optional[Match]
    for item in structures:
        summary = item.get("summary", item) if isinstance(item, dict) else {}
        match = _afdb_match_from_sequence_summary(record, summary)
        if match is None:
            continue
        if best is None or _is_better_match(match, best):
            best = match
    return best


def _afdb_match_from_sequence_summary(record: SequenceRecord, summary: Dict[str, Any]) -> Optional[Match]:
    model_id = summary.get("model_identifier")
    if not model_id:
        return None
    seq_id = _optional_float(summary.get("sequence_identity"))
    coverage = _optional_float(summary.get("coverage"))
    if seq_id is None or coverage is None:
        return None
    seq_id = _as_fraction(seq_id)
    coverage = _as_fraction(coverage)
    model_url = _afdb_model_url_from_sequence_summary(summary)
    return Match(
        protein_id=record.protein_id,
        match="AFDB:{}".format(model_id),
        evalue=None,
        seq_id=seq_id,
        query_coverage=coverage,
        target_coverage=coverage,
        query_length=len(record.sequence),
        target_length=len(record.sequence),
        align_length=max(1, int(round(len(record.sequence) * coverage))),
        plddt=_optional_float(summary.get("confidence_avg_local_score")),
        pdb_path=model_url,
    )


def _afdb_model_url_from_sequence_summary(summary: Dict[str, Any]) -> Optional[str]:
    model_url = summary.get("model_url")
    if not model_url:
        return None
    return model_url


def _run_diamond(
    query: Path,
    db: Path,
    output: Path,
    threads: int,
    diamond_executable: Optional[str] = None,
) -> None:
    diamond = diamond_executable or shutil.which("diamond")
    if diamond is None:
        raise FileNotFoundError(
            "diamond executable not found; set decision_tree.diamond_executable, "
            "set decision_tree.pdb_hits_tsv, or install DIAMOND on PATH"
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    subprocess.check_call(
        [
            diamond,
            "blastp",
            "--db",
            str(db),
            "--query",
            str(query),
            "-o",
            str(output),
            "--threads",
            str(threads),
            "--fast",
            "--outfmt",
            "6",
            "qseqid",
            "sseqid",
            "evalue",
            "nident",
            "qlen",
            "slen",
            "length",
            "--top",
            "1",
        ]
    )


def _decision_match(match: Match) -> str:
    if match.pdb_path and _is_url(match.pdb_path):
        return "{}|model_url={}".format(match.match, quote(match.pdb_path, safe=":/?&=.-_~"))
    return match.match


def _read_hits_tsv(path: Path, source: str) -> Dict[str, Match]:
    _require_file(path, "{} hit report".format(source))
    best = {}  # type: Dict[str, Match]
    with path.open(encoding="utf-8") as handle:
        first = handle.readline()
        if not first:
            return {}
        handle.seek(0)
        has_header = any(name in first.lower().split("\t") for name in ("protein_id", "qseqid", "query"))
        if has_header:
            reader = csv.DictReader(handle, delimiter="\t")
            rows = (_match_from_mapping(row, source) for row in reader)
        else:
            rows = (_match_from_columns(line.rstrip("\n").split("\t"), source) for line in handle if line.strip())
        for match in rows:
            current = best.get(match.protein_id)
            if current is None or _is_better_match(match, current):
                best[match.protein_id] = match
    return best


def _match_from_mapping(row: Dict[str, str], source: str) -> Match:
    protein_id = _protein_id(row.get("protein_id") or row.get("query") or row.get("qseqid") or "")
    match = row.get("match") or row.get("target") or row.get("sseqid") or ""
    evalue = _optional_float(row.get("evalue") or row.get("E"))
    nident = _optional_int(row.get("nident"))
    qlen = _required_int(row.get("qlen") or row.get("query_length"))
    slen = _required_int(row.get("slen") or row.get("target_length"))
    align_length = _optional_int(row.get("length") or row.get("align_length")) or nident or min(qlen, slen)
    seq_id = _optional_float(row.get("seq_id") or row.get("identity"))
    if seq_id is None:
        if nident is None:
            raise ValueError("Hit report row needs seq_id/identity or nident: {}".format(row))
        seq_id = float(nident) / float(align_length)
    seq_id = _as_fraction(seq_id)
    query_coverage = _optional_float(row.get("query_coverage") or row.get("qcov"))
    target_coverage = _optional_float(row.get("target_coverage") or row.get("tcov"))
    if query_coverage is None:
        query_coverage = float(align_length) / float(qlen)
    if target_coverage is None:
        target_coverage = float(align_length) / float(slen)
    query_coverage = _as_fraction(query_coverage)
    target_coverage = _as_fraction(target_coverage)
    return Match(
        protein_id=protein_id,
        match=_normalize_match(match, source),
        evalue=evalue,
        seq_id=seq_id,
        query_coverage=query_coverage,
        target_coverage=target_coverage,
        query_length=qlen,
        target_length=slen,
        align_length=align_length,
        plddt=_optional_float(row.get("plddt") or row.get("mean_plddt")),
        pdb_path=row.get("pdb_path") or None,
    )


def _match_from_columns(parts: List[str], source: str) -> Match:
    if len(parts) < 6:
        raise ValueError("Expected at least 6 tab-delimited columns in hit report")
    qseqid, sseqid, evalue, nident, qlen, slen = parts[:6]
    align_length = int(parts[6]) if len(parts) > 6 and parts[6] else int(nident)
    qlen_i = int(qlen)
    slen_i = int(slen)
    nident_i = int(nident)
    return Match(
        protein_id=_protein_id(qseqid),
        match=_normalize_match(sseqid, source),
        evalue=float(evalue),
        seq_id=float(nident_i) / float(align_length),
        query_coverage=float(align_length) / float(qlen_i),
        target_coverage=float(align_length) / float(slen_i),
        query_length=qlen_i,
        target_length=slen_i,
        align_length=align_length,
    )


def mean_plddt_from_pdb(path: Path) -> Optional[float]:
    total = 0.0
    count = 0
    seen = set()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            residue = (line[21:22], line[22:26].strip(), line[26:27].strip())
            if residue in seen:
                continue
            seen.add(residue)
            try:
                total += float(line[60:66].strip())
            except ValueError:
                continue
            count += 1
    if count == 0:
        return None
    return total / float(count)


def _has_good_plddt(match: Match, threshold: float) -> bool:
    return match.plddt is not None and match.plddt >= threshold


def _score(value: float) -> str:
    return "{:.6g}".format(value)


def _is_better_match(candidate: Match, current: Match) -> bool:
    if candidate.evalue is not None and current.evalue is not None and candidate.evalue != current.evalue:
        return candidate.evalue < current.evalue
    candidate_plddt = candidate.plddt if candidate.plddt is not None else -1.0
    current_plddt = current.plddt if current.plddt is not None else -1.0
    return (candidate.seq_id, candidate.query_coverage, candidate_plddt) > (
        current.seq_id,
        current.query_coverage,
        current_plddt,
    )


def _normalize_match(match: str, source: str) -> str:
    if source == "AFDB" and match:
        if match.startswith("AFDB:"):
            return "AFDB:{}".format(_strip_afdb_model_suffix(match.split(":", 1)[1]))
        return "AFDB:{}".format(_strip_afdb_model_suffix(match))
    return match


def _strip_afdb_model_suffix(match: str) -> str:
    for suffix in (".pdb", ".cif", ".bcif"):
        if match.endswith(suffix):
            match = match[: -len(suffix)]
            break
    marker = "-model_v"
    if marker in match:
        return match.split(marker, 1)[0]
    return match


def _protein_id(value: str) -> str:
    if value.startswith(">"):
        return protein_id_from_header(value)
    parts = value.split("|")
    if len(parts) >= 2:
        return parts[1]
    return value


def _pdb_dmnd_path(decision_cfg: Dict[str, Any]) -> Optional[Path]:
    if decision_cfg.get("pdb_dmnd"):
        return Path(decision_cfg["pdb_dmnd"]).expanduser()
    return None


def _is_url(value: str) -> bool:
    return value.startswith("http://") or value.startswith("https://")


def _clean_sequence(sequence: str) -> str:
    return "".join(sequence.split()).upper()


def _pdb_search_plan(config: Dict[str, Any]) -> Dict[str, Any]:
    decision_cfg = config.get("decision_tree", {})
    return {
        "hits_tsv": decision_cfg.get("pdb_hits_tsv"),
        "pdb_dmnd": str(_pdb_dmnd_path(decision_cfg)) if _pdb_dmnd_path(decision_cfg) else None,
    }


def _afdb_search_plan(config: Dict[str, Any]) -> Dict[str, Any]:
    decision_cfg = config.get("decision_tree", {})
    return {
        "hits_tsv": decision_cfg.get("afdb_hits_tsv"),
        "use_afdb_api": decision_cfg.get("use_afdb_api", True),
        "sequence_search": decision_cfg.get("afdb_sequence_search", True),
    }


def _optional_float(value: Any) -> Optional[float]:
    if value in (None, ""):
        return None
    return float(value)


def _as_fraction(value: float) -> float:
    if value > 1.0:
        return value / 100.0
    return value


def _first_float(mapping: Dict[str, Any], keys: Iterable[str]) -> Optional[float]:
    for key in keys:
        value = mapping.get(key)
        if value not in (None, ""):
            try:
                return float(value)
            except (TypeError, ValueError):
                pass
    return None


def _first_present(mapping: Dict[str, Any], keys: Iterable[str]) -> Optional[str]:
    for key in keys:
        value = mapping.get(key)
        if value not in (None, ""):
            return str(value)
    return None


def _optional_int(value: Any) -> Optional[int]:
    if value in (None, ""):
        return None
    return int(float(value))


def _required_int(value: Any) -> int:
    if value in (None, ""):
        raise ValueError("Missing required integer field in hit report")
    return int(float(value))


def _looks_like_uniprot_accession(value: str) -> bool:
    return 6 <= len(value) <= 10 and value.isalnum()


def _write_metadata(run: DecisionTreeRun) -> None:
    metadata = {
        "stage": "decision_tree",
        "mode": "native",
        "run_dir": str(run.run_dir),
        "work_dir": str(run.work_dir),
        "pipeline_files": str(run.pipeline_files),
        "input_fasta": str(run.input_fasta),
        "decisions": str(run.decisions),
        "command": run.command,
        "command_display": format_command(run.command, run.work_dir),
        "pdb_search": _pdb_search_plan(run.config),
        "afdb_search": _afdb_search_plan(run.config),
    }
    (run.run_dir / "decision_tree_run.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )


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
        raise FileNotFoundError("Missing {}: {}".format(label, path))
