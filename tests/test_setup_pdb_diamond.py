import tempfile
import unittest
from pathlib import Path

from scripts.setup_pdb_diamond import _read_protein_entries, _write_protein_fasta


class SetupPdbDiamondTests(unittest.TestCase):
    def test_filters_pdb_seqres_to_protein_records(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            entry_types = root / "pdb_entry_type.txt"
            entry_types.write_text(
                "1abc prot diffraction\n"
                "2xyz nuc diffraction\n",
                encoding="utf-8",
            )
            seqres = root / "pdb_seqres.txt"
            seqres.write_text(
                ">1abc_A mol:protein length:3\n"
                "MKT\n"
                ">2xyz_A mol:na length:3\n"
                "AAA\n",
                encoding="utf-8",
            )
            output = root / "pdb.fasta"

            proteins = _read_protein_entries(entry_types)
            kept = _write_protein_fasta(seqres, proteins, output)

            self.assertEqual(kept, 1)
            self.assertEqual(
                output.read_text(encoding="utf-8"),
                ">1abc_A mol:protein length:3\nMKT\n",
            )


if __name__ == "__main__":
    unittest.main()
