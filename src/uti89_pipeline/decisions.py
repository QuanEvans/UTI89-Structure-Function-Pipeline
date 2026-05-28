"""Decision-file parsing."""

from pathlib import Path
from typing import Dict, Iterator, NamedTuple, Optional, Tuple


class Decision(NamedTuple):
    protein_id: str
    method: str
    match: str
    score: str
    reason: str


def read_decisions(path: Path) -> Iterator[Decision]:
    """Yield decision rows from the decision-tree TSV output."""
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                raise ValueError(
                    "Expected at least 5 tab-delimited columns in {} line {}".format(
                        path, line_number
                    )
                )
            yield Decision(
                protein_id=parts[0],
                method=parts[1],
                match=parts[2],
                score=parts[3],
                reason=parts[4],
            )


def read_decision_map(path: Path) -> Dict[str, Decision]:
    """Return decisions keyed by protein ID."""
    return {decision.protein_id: decision for decision in read_decisions(path)}


def afdb_model_url(match: str) -> Optional[str]:
    """Return an AFDB model URL carried in a decision match field, if present."""
    _base_match, metadata = split_match_metadata(match)
    return metadata.get("model_url") or metadata.get("pdb_url")


def split_match_metadata(match: str) -> Tuple[str, Dict[str, str]]:
    """Split a match field into its primary match token and optional key/value metadata."""
    if not match:
        return "", {}
    parts = match.split("|")
    metadata = {}  # type: Dict[str, str]
    for item in parts[1:]:
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        metadata[key] = value
    return parts[0], metadata
