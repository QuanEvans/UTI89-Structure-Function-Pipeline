#!/usr/bin/env python3
"""Download PDB SEQRES proteins and build the DIAMOND database used by the decision tree."""

import argparse
import subprocess
import urllib.request
from pathlib import Path


SEQRES_URLS = [
    "https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt",
    "https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt",
]
ENTRY_TYPE_URLS = [
    "https://files.wwpdb.org/pub/pdb/derived_data/pdb_entry_type.txt",
    "https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_entry_type.txt",
]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out-dir",
        required=True,
        help="Directory where pdb_seqres.txt, pdb.fasta, and pdb.dmnd will be written.",
    )
    parser.add_argument(
        "--diamond",
        default="diamond",
        help="DIAMOND executable to use for makedb.",
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Reuse existing pdb_seqres.txt and pdb_entry_type.txt in --out-dir.",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=60.0,
        help="Seconds to wait for each wwPDB download URL before trying the next mirror.",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    seqres = out_dir / "pdb_seqres.txt"
    entry_types = out_dir / "pdb_entry_type.txt"
    protein_fasta = out_dir / "pdb.fasta"
    db_prefix = out_dir / "pdb"

    if not args.skip_download:
        _download_first_available(SEQRES_URLS, seqres, args.timeout)
        _download_first_available(ENTRY_TYPE_URLS, entry_types, args.timeout)

    if not seqres.is_file():
        raise FileNotFoundError("Missing PDB SEQRES FASTA: {}".format(seqres))
    if not entry_types.is_file():
        raise FileNotFoundError("Missing PDB entry type file: {}".format(entry_types))

    protein_entries = _read_protein_entries(entry_types)
    kept = _write_protein_fasta(seqres, protein_entries, protein_fasta)
    if kept == 0:
        raise RuntimeError("No protein records were written from {}".format(seqres))

    subprocess.check_call(
        [
            args.diamond,
            "makedb",
            "--in",
            str(protein_fasta),
            "-d",
            str(db_prefix),
        ]
    )
    print("Wrote {}".format(protein_fasta))
    print("Wrote {}.dmnd".format(db_prefix))
    print("Set decision_tree.pdb_dmnd: {}.dmnd".format(db_prefix))
    return 0


def _download_first_available(urls: list, destination: Path, timeout: float) -> None:
    errors = []
    for url in urls:
        try:
            _download(url, destination, timeout)
            return
        except Exception as exc:
            errors.append("{}: {}".format(url, exc))
    raise RuntimeError(
        "Could not download {} from any configured URL:\n{}".format(
            destination.name, "\n".join(errors)
        )
    )


def _download(url: str, destination: Path, timeout: float) -> None:
    tmp = destination.with_suffix(destination.suffix + ".tmp")
    print("Downloading {} -> {}".format(url, destination))
    with urllib.request.urlopen(url, timeout=timeout) as response, tmp.open("wb") as handle:
        handle.write(response.read())
    tmp.replace(destination)


def _read_protein_entries(path: Path) -> set:
    proteins = set()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2 and "prot" in parts[1].lower():
                proteins.add(parts[0].lower())
    return proteins


def _write_protein_fasta(seqres: Path, protein_entries: set, output: Path) -> int:
    kept = 0
    write_record = False
    with seqres.open(encoding="utf-8") as source, output.open("w", encoding="utf-8") as dest:
        for line in source:
            if line.startswith(">"):
                entry = line[1:5].lower()
                header_lower = line.lower()
                write_record = entry in protein_entries or "mol:protein" in header_lower
                if write_record:
                    kept += 1
                    dest.write(line)
                continue
            if write_record:
                dest.write(line)
    return kept


if __name__ == "__main__":
    raise SystemExit(main())
