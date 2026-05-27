import unittest

from scripts.run_pipeline import _apply_work_dir_override, _uses_user_models


class RunPipelineTests(unittest.TestCase):
    def test_work_dir_override_clears_derived_stage_paths(self) -> None:
        config = {
            "run": {"work_dir": "/old"},
            "structure_prediction": {
                "output_dir": "/old/structures",
                "decision_file": "/old/decisions.txt",
                "afdb_pdb_dir": "/old/pdb",
                "enabled_steps": ["AFDB"],
            },
            "function_prediction": {
                "structure_dir": "/old/structures",
                "output_dir": "/old/functions",
                "decision_file": "/old/decisions.txt",
                "starfunc_sif": "/example/StarFunc.sif",
            },
        }

        _apply_work_dir_override(config, "/new")

        self.assertEqual(config["run"]["work_dir"], "/new")
        self.assertEqual(config["structure_prediction"]["enabled_steps"], ["AFDB"])
        self.assertEqual(config["function_prediction"]["starfunc_sif"], "/example/StarFunc.sif")
        for key in ("output_dir", "decision_file", "afdb_pdb_dir"):
            self.assertNotIn(key, config["structure_prediction"])
        for key in ("structure_dir", "output_dir", "decision_file"):
            self.assertNotIn(key, config["function_prediction"])

    def test_uses_user_models(self) -> None:
        self.assertTrue(
            _uses_user_models({"structure_prediction": {"user_models": {"model_dir": "/models"}}})
        )
        self.assertFalse(
            _uses_user_models(
                {"structure_prediction": {"user_models": {"enabled": False, "model_dir": "/models"}}}
            )
        )


if __name__ == "__main__":
    unittest.main()
