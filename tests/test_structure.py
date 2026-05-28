import tempfile
import unittest
from pathlib import Path

from uti89_pipeline.structure import prepare_structure_predictions


class StructureTests(unittest.TestCase):
    def test_prepare_structure_predictions(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            repo = Path.cwd()
            config = {
                "run": {
                    "work_dir": str(root / "run"),
                    "input_fasta": str(repo / "tests" / "data" / "mini.fasta"),
                },
                "slurm": {"partition": "compute", "account": "project_account"},
                "structure_prediction": {
                    "decision_file": str(repo / "tests" / "data" / "decisions.tsv"),
                    "output_dir": str(root / "structures"),
                    "enabled_steps": ["GETSEQS", "AFDB", "DITASSER"],
                    "tools": {
                        "seq2fun_sif": "/example/seq2fun.sif",
                        "ditasser_pkgdir": "/example/D-I-TASSER-2.0",
                        "itlib_dir": "/example/ITLIB",
                        "scratch_bind": "/example/scratch",
                    },
                },
            }

            result = prepare_structure_predictions(config)

            self.assertEqual(len(result.prepared), 2)
            self.assertTrue((root / "structures" / "P00001" / "seq.fasta").is_file())
            self.assertTrue((root / "structures" / "P00001" / "fetch_afdb.sbatch").is_file())
            self.assertTrue((root / "structures" / "P00002" / "run_ditasser.sbatch").is_file())

    def test_prepare_user_model_structures(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            repo = Path.cwd()
            model_dir = root / "models"
            model_dir.mkdir()
            (model_dir / "P00001.pdb").write_text("HEADER P00001\nEND\n", encoding="utf-8")
            (model_dir / "P00002.pdb").write_text("HEADER P00002\nEND\n", encoding="utf-8")
            config = {
                "run": {
                    "work_dir": str(root / "run"),
                    "input_fasta": str(repo / "tests" / "data" / "mini.fasta"),
                },
                "structure_prediction": {
                    "output_dir": str(root / "structures"),
                    "user_models": {
                        "model_dir": str(model_dir),
                        "pattern": "{protein_id}.pdb",
                        "link_models": False,
                    },
                },
            }

            result = prepare_structure_predictions(config)

            self.assertEqual(len(result.prepared), 2)
            self.assertFalse(result.skipped)
            self.assertTrue((root / "structures" / "P00001" / "model1.pdb").is_file())
            self.assertTrue((root / "structures" / "P00002" / "seq.fasta").is_file())

    def test_prepare_modeler_fetches_afdb_template_url(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            repo = Path.cwd()
            decisions = root / "decisions.tsv"
            decisions.write_text(
                "P00001\tMODELLER\t"
                "AFDB:AF-P00001-F1|model_url=https://example.org/AF-P00001-F1-model_v6.cif"
                "\t0.96\tafdb_template_high_confidence\n",
                encoding="utf-8",
            )
            config = {
                "run": {
                    "work_dir": str(root / "run"),
                    "input_fasta": str(repo / "tests" / "data" / "mini.fasta"),
                },
                "structure_prediction": {
                    "decision_file": str(decisions),
                    "output_dir": str(root / "structures"),
                    "enabled_steps": ["GETSEQS", "MODELLER"],
                    "tools": {
                        "seq2fun_sif": "/example/seq2fun.sif",
                        "ditasser_pkgdir": "/example/D-I-TASSER-2.0",
                        "itlib_dir": "/example/ITLIB",
                        "scratch_bind": "/example/scratch",
                    },
                },
            }

            result = prepare_structure_predictions(config)

            self.assertEqual(len(result.prepared), 1)
            script = root / "structures" / "P00001" / "run_modeler.sh"
            text = script.read_text(encoding="utf-8")
            self.assertIn("urllib.request.urlretrieve", text)
            self.assertIn("https://example.org/AF-P00001-F1-model_v6.cif", text)
            self.assertIn("import gemmi", text)
            self.assertIn("gemmi.read_structure", text)
            self.assertIn('suffix in (".cif", ".bcif")', text)
            self.assertIn("template.pdb", text)

    def test_prepare_afdb_fetch_script_uses_api_model_url(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            repo = Path.cwd()
            decisions = root / "decisions.tsv"
            decisions.write_text(
                "P00001\tAFDB\t"
                "AFDB:AF-P00001-F1|model_url=https://example.org/AF-P00001-F1-model_v6.cif"
                "\t1\tafdb_exact_high_confidence\n",
                encoding="utf-8",
            )
            config = {
                "run": {
                    "work_dir": str(root / "run"),
                    "input_fasta": str(repo / "tests" / "data" / "mini.fasta"),
                },
                "structure_prediction": {
                    "decision_file": str(decisions),
                    "output_dir": str(root / "structures"),
                    "enabled_steps": ["GETSEQS", "AFDB"],
                },
            }

            result = prepare_structure_predictions(config)

            self.assertEqual(len(result.prepared), 1)
            script = root / "structures" / "P00001" / "fetch_afdb.sh"
            text = script.read_text(encoding="utf-8")
            self.assertIn("cd {}".format(root / "structures" / "P00001"), text)
            self.assertIn("urllib.request.urlretrieve", text)
            self.assertIn("https://example.org/AF-P00001-F1-model_v6.cif", text)
            self.assertIn("import gemmi", text)
            self.assertIn("gemmi.read_structure", text)
            self.assertIn('suffix in (".cif", ".bcif")', text)
            self.assertIn("model1.pdb", text)
