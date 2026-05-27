"""FASTA helpers."""

from pathlib import Path
from typing import Iterator, Tuple


def read_fasta_records(path: Path) -> Iterator[Tuple[str, str]]:
    """Yield `(header, sequence)` records from a FASTA file."""
    header = None
    seq_lines = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line
                seq_lines = []
                continue
            if header is None:
                raise ValueError("Found sequence data before first FASTA header")
            seq_lines.append(line)
    if header is not None:
        yield header, "".join(seq_lines)


def protein_id_from_header(header: str) -> str:
    """Extract the protein ID used by the decision-tree output."""
    clean = header[1:] if header.startswith(">") else header
    parts = clean.split("|")
    if len(parts) >= 2 and parts[1]:
        return parts[1]
    return clean.split()[0]


def write_single_record(path: Path, protein_id: str, sequence: str) -> None:
    """Write one FASTA record with a normalized protein ID header."""
    path.write_text(">{}\n{}\n".format(protein_id, sequence), encoding="utf-8")

