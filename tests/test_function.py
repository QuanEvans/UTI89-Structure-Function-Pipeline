import tempfile
import unittest
from pathlib import Path

from uti89_pipeline.function import prepare_function_predictions


class FunctionTests(unittest.TestCase):
    def test_prepare_function_predictions(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            repo = Path.cwd()
            structure_dir = root / "structures" / "P00001"
            structure_dir.mkdir(parents=True)
            (structure_dir / "model1.pdb").write_text("HEADER TEST\nEND\n", encoding="utf-8")
            (structure_dir / "seq.fasta").write_text(">P00001\nMKT\n", encoding="utf-8")
            config = {
                "run": {"work_dir": str(root / "run"), "input_fasta": "unused.fasta"},
                "execution": {"backend": "slurm"},
                "slurm": {"partition": "compute", "account": "project_account"},
                "structure_prediction": {"output_dir": str(root / "structures")},
                "function_prediction": {
                    "decision_file": str(repo / "tests" / "data" / "decisions.tsv"),
                    "output_dir": str(root / "functions"),
                    "starfunc_sif": "/example/StarFunc.sif",
                    "starfunc_database": "/example/database",
                },
            }

            result = prepare_function_predictions(config)

            self.assertEqual(len(result.prepared), 1)
            self.assertTrue((root / "functions" / "P00001" / "input.pdb").exists())
            self.assertTrue((root / "functions" / "P00001" / "run_starfunc.sbatch").is_file())
            self.assertIn("P00002: missing model1.pdb", result.skipped)

    def test_prepare_function_predictions_without_decisions(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            structure_dir = root / "structures" / "P00001"
            structure_dir.mkdir(parents=True)
            (structure_dir / "model1.pdb").write_text("HEADER TEST\nEND\n", encoding="utf-8")
            (structure_dir / "seq.fasta").write_text(">P00001\nMKT\n", encoding="utf-8")
            config = {
                "run": {"work_dir": str(root / "run"), "input_fasta": "unused.fasta"},
                "structure_prediction": {"output_dir": str(root / "structures")},
                "function_prediction": {
                    "output_dir": str(root / "functions"),
                    "starfunc_sif": "/example/StarFunc.sif",
                    "starfunc_database": "/example/database",
                },
            }

            result = prepare_function_predictions(config)

            self.assertEqual(len(result.prepared), 1)
            self.assertFalse(result.skipped)
            script = root / "functions" / "P00001" / "run_starfunc.sh"
            self.assertTrue(script.is_file())
            self.assertNotIn("#SBATCH", script.read_text(encoding="utf-8"))
