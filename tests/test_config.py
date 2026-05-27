import unittest
from pathlib import Path

from uti89_pipeline.config import load_config


class ConfigTests(unittest.TestCase):
    def test_load_simple_yaml_example(self) -> None:
        config = load_config(Path("config/template.yaml"))

        self.assertEqual(config["execution"]["backend"], "local")
        self.assertEqual(config["decision_tree"]["af_min_seq_id"], 0.9)
        self.assertEqual(
            config["structure_prediction"]["enabled_steps"],
            ["GETSEQS", "AFDB", "MODELLER", "LOMETS", "DITASSER"],
        )


if __name__ == "__main__":
    unittest.main()
