import unittest

from uti89_pipeline.execution import backend, module_load_line, script_suffix


class ExecutionTests(unittest.TestCase):
    def test_defaults_to_local_without_slurm_section(self) -> None:
        config = {"execution": {"local_cores": 4}}

        self.assertEqual(backend(config), "local")
        self.assertEqual(script_suffix(config), ".sh")
        self.assertEqual(module_load_line(config, "singularity"), "")

    def test_slurm_backend_is_explicit_or_legacy_slurm_section(self) -> None:
        config = {
            "execution": {"backend": "slurm"},
            "slurm": {"partition": "compute", "account": "project_account"},
        }

        self.assertEqual(backend(config), "slurm")
        self.assertEqual(script_suffix(config), ".sbatch")
        self.assertIn("module load singularity", module_load_line(config, "singularity"))


if __name__ == "__main__":
    unittest.main()
