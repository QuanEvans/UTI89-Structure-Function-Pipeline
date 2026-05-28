import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scripts.setup_pdb_diamond import _download_first_available, _read_protein_entries, _write_protein_fasta


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

    def test_download_uses_timeout(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            output = Path(tempdir) / "pdb_seqres.txt"
            with patch("scripts.setup_pdb_diamond.urllib.request.urlopen") as urlopen:
                urlopen.return_value.__enter__.return_value.read.return_value = b">1abc_A\nMKT\n"

                _download_first_available(["https://example.org/pdb_seqres.txt"], output, 12.5)

            urlopen.assert_called_once_with("https://example.org/pdb_seqres.txt", timeout=12.5)
            self.assertEqual(output.read_bytes(), b">1abc_A\nMKT\n")


if __name__ == "__main__":
    unittest.main()
