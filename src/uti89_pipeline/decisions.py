"""Decision-file parsing."""

from pathlib import Path
from typing import Dict, Iterator, NamedTuple, Optional


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


def afdb_model_filename(match: str) -> Optional[str]:
    """Return the AFDB model filename implied by a decision match field."""
    if not match or not match.startswith("AFDB:"):
        return None
    af_name = match.split(":", 1)[1]
    return "{}-model_v4.pdb".format(af_name)

