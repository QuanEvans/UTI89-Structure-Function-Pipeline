import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from uti89_pipeline.decision_tree import (
    Match,
    choose_structure_method,
    prepare_decision_tree_run,
    run_decision_tree,
)


class DecisionTreeTests(unittest.TestCase):
    def test_native_decision_tree_uses_afdb_hit_report_url(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            repo = Path.cwd()
            afdb_hits = root / "afdb_hits.tsv"
            afdb_hits.write_text(
                "protein_id\tmatch\tseq_id\tquery_coverage\ttarget_coverage\tqlen\tslen\tplddt\tpdb_path\n"
                "P00001\tAFDB:AF-P00001-F1\t1.0\t1.0\t1.0\t26\t26\t91\t"
                "https://alphafold.ebi.ac.uk/files/AF-P00001-F1-model_v6.pdb\n",
                encoding="utf-8",
            )
            config = {
                "run": {
                    "work_dir": str(root / "run"),
                    "input_fasta": str(repo / "tests" / "data" / "mini.fasta"),
                },
                "decision_tree": {
                    "afdb_hits_tsv": str(afdb_hits),
                    "use_afdb_api": False,
                },
            }

            run = prepare_decision_tree_run(config)
            self.assertEqual(run_decision_tree(run), 0)

            decisions = run.decisions.read_text(encoding="utf-8").splitlines()
            self.assertEqual(
                decisions[0],
                "P00001\tAFDB\t"
                "AFDB:AF-P00001-F1|model_url=https://alphafold.ebi.ac.uk/files/AF-P00001-F1-model_v6.pdb"
                "\t1\tafdb_exact_high_confidence",
            )
            self.assertEqual(
                decisions[1],
                "P00002\tDITASSER\t\t0\tno_high_confidence_template",
            )

    def test_native_decision_tree_uses_hit_reports(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            repo = Path.cwd()
            afdb_hits = root / "afdb_hits.tsv"
            afdb_hits.write_text(
                "protein_id\tmatch\tseq_id\tquery_coverage\ttarget_coverage\tqlen\tslen\tplddt\n"
                "P00001\tAFDB:AF-P00001-F1\t0.96\t1.0\t1.0\t26\t26\t88\n",
                encoding="utf-8",
            )
            pdb_hits = root / "pdb_hits.tsv"
            pdb_hits.write_text(
                "protein_id\tmatch\tseq_id\tquery_coverage\ttarget_coverage\tqlen\tslen\n"
                "P00002\t1ABC_A\t0.98\t0.95\t0.95\t20\t20\n",
                encoding="utf-8",
            )
            config = {
                "run": {
                    "work_dir": str(root / "run"),
                    "input_fasta": str(repo / "tests" / "data" / "mini.fasta"),
                },
                "decision_tree": {
                    "afdb_hits_tsv": str(afdb_hits),
                    "pdb_hits_tsv": str(pdb_hits),
                    "use_afdb_api": False,
                },
            }

            run = prepare_decision_tree_run(config)
            self.assertEqual(run_decision_tree(run), 0)

            decisions = run.decisions.read_text(encoding="utf-8").splitlines()
            self.assertEqual(
                decisions[0],
                "P00001\tMODELLER\tAFDB:AF-P00001-F1\t0.96\tafdb_template_high_confidence",
            )
            self.assertEqual(
                decisions[1],
                "P00002\tLOMETS\t1ABC_A\t0.98\tpdb_high_identity",
            )

    def test_native_decision_tree_uses_afdb_api_for_exact_accession(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            fasta = root / "input.fasta"
            fasta.write_text(">P07550\nMTESTSEQUENCE\n", encoding="utf-8")
            config = {
                "run": {
                    "work_dir": str(root / "run"),
                    "input_fasta": str(fasta),
                },
                "decision_tree": {"use_afdb_api": True},
            }

            payload = [
                {
                    "entryId": "AF-P07550-F1",
                    "globalMetricValue": 88.5,
                    "sequence": "MTESTSEQUENCE",
                    "cifUrl": "https://alphafold.ebi.ac.uk/files/AF-P07550-F1-model_v6.cif",
                    "pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-P07550-F1-model_v6.pdb",
                }
            ]
            with patch("uti89_pipeline.decision_tree.urllib.request.urlopen") as urlopen:
                urlopen.return_value = _FakeResponse(payload)
                run = prepare_decision_tree_run(config)
                self.assertEqual(run_decision_tree(run), 0)

            decisions = run.decisions.read_text(encoding="utf-8").splitlines()
            self.assertEqual(
                decisions[0],
                "P07550\tAFDB\t"
                "AFDB:AF-P07550-F1|model_url=https://alphafold.ebi.ac.uk/files/AF-P07550-F1-model_v6.cif"
                "\t1\tafdb_exact_high_confidence",
            )

    def test_native_decision_tree_rejects_afdb_api_sequence_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            fasta = root / "input.fasta"
            fasta.write_text(">P07550\nMQUERY\n", encoding="utf-8")
            config = {
                "run": {
                    "work_dir": str(root / "run"),
                    "input_fasta": str(fasta),
                },
                "decision_tree": {"use_afdb_api": True},
            }

            payload = [{"entryId": "AF-P07550-F1", "globalMetricValue": 88.5, "sequence": "MDIFFERENT"}]
            with patch("uti89_pipeline.decision_tree.urllib.request.urlopen") as urlopen:
                urlopen.return_value = _FakeResponse(payload)
                run = prepare_decision_tree_run(config)
                self.assertEqual(run_decision_tree(run), 0)

            decisions = run.decisions.read_text(encoding="utf-8").splitlines()
            self.assertEqual(decisions[0], "P07550\tDITASSER\t\t0\tno_high_confidence_template")

    def test_native_decision_tree_uses_afdb_sequence_search_without_uniprot_id(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            sequence = "MTESTSEQAAMTESTSEQAA"
            fasta = root / "input.fasta"
            fasta.write_text(">my_protein\n{}\n".format(sequence), encoding="utf-8")
            config = {
                "run": {
                    "work_dir": str(root / "run"),
                    "input_fasta": str(fasta),
                },
                "decision_tree": {"use_afdb_api": True, "afdb_sequence_search": True},
            }

            payload = {
                "entry": {"sequence": sequence, "checksum": "abc", "checksum_type": "MD5"},
                "structures": [
                    {
                        "summary": {
                            "model_identifier": "AF-Q11111-F1",
                            "model_url": "https://alphafold.ebi.ac.uk/files/AF-Q11111-F1-model_v6.cif",
                            "sequence_identity": 1.0,
                            "coverage": 1.0,
                            "confidence_avg_local_score": 86.0,
                        }
                    }
                ],
            }
            with patch("uti89_pipeline.decision_tree.urllib.request.urlopen") as urlopen:
                urlopen.return_value = _FakeResponse(payload)
                run = prepare_decision_tree_run(config)
                self.assertEqual(run_decision_tree(run), 0)

            called_url = urlopen.call_args[0][0]
            self.assertIn("/api/sequence/summary?", called_url)
            self.assertIn("type=sequence", called_url)
            decisions = run.decisions.read_text(encoding="utf-8").splitlines()
            self.assertEqual(
                decisions[0],
                "my_protein\tAFDB\t"
                "AFDB:AF-Q11111-F1|model_url=https://alphafold.ebi.ac.uk/files/AF-Q11111-F1-model_v6.cif"
                "\t1\tafdb_exact_high_confidence",
            )

    def test_native_decision_tree_can_run_pdb_diamond_search(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            repo = Path.cwd()
            bin_dir = root / "bin"
            bin_dir.mkdir()
            diamond = bin_dir / "diamond"
            diamond.write_text(
                "#!/usr/bin/env python3\n"
                "import sys\n"
                "out = sys.argv[sys.argv.index('-o') + 1]\n"
                "open(out, 'w', encoding='utf-8').write("
                "'P00002\\t1ABC_A\\t1e-50\\t20\\t20\\t20\\t20\\n')\n",
                encoding="utf-8",
            )
            diamond.chmod(0o755)
            pdb_dmnd = root / "pdb.dmnd"
            pdb_dmnd.write_text("placeholder\n", encoding="utf-8")
            config = {
                "run": {
                    "work_dir": str(root / "run"),
                    "input_fasta": str(repo / "tests" / "data" / "mini.fasta"),
                },
                "decision_tree": {
                    "pdb_dmnd": str(pdb_dmnd),
                    "diamond_executable": str(diamond),
                    "use_afdb_api": False,
                },
            }

            run = prepare_decision_tree_run(config)
            self.assertEqual(run_decision_tree(run), 0)

            decisions = run.decisions.read_text(encoding="utf-8").splitlines()
            self.assertEqual(decisions[1], "P00002\tLOMETS\t1ABC_A\t1\tpdb_high_identity")

    def test_native_decision_tree_replays_representative_original_decisions(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            fasta = root / "input.fasta"
            fasta.write_text(
                "\n".join(
                    [
                        ">Q1R393",
                        "M" * 60,
                        ">A0A080WMA9",
                        "M" * 60,
                        ">A0A075B5P5",
                        "M" * 60,
                        ">P03690",
                        "M" * 60,
                        ">Q1R415",
                        "M" * 60,
                        ">A0A7S9SVH6",
                        "M" * 60,
                        ">A0A7S9SVS1",
                        "M" * 60,
                        ">H9L4R6",
                        "M" * 60,
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            afdb_hits = root / "afdb_hits.tsv"
            afdb_hits.write_text(
                "protein_id\tmatch\tseq_id\tquery_coverage\ttarget_coverage\tqlen\tslen\tplddt\tpdb_path\n"
                "Q1R393\tAFDB:AF-A0A0H2VDW3-F1\t1.0\t1.0\t1.0\t60\t60\t90\t"
                "https://alphafold.ebi.ac.uk/files/AF-A0A0H2VDW3-F1-model_v6.pdb\n"
                "A0A080WMA9\tAFDB:AF-A0A080WMA9-F1\t1.0\t1.0\t1.0\t60\t60\t90\t"
                "https://alphafold.ebi.ac.uk/files/AF-A0A080WMA9-F1-model_v6.pdb\n"
                "A0A075B5P5\tAFDB:AF-A0A4U9FKB1-F1\t0.996969696969697\t1.0\t1.0\t60\t60\t88\t"
                "https://alphafold.ebi.ac.uk/files/AF-A0A4U9FKB1-F1-model_v6.pdb\n"
                "P03690\tAFDB:AF-A0A2M9X344-F1\t0.9806896551724138\t1.0\t1.0\t60\t60\t88\t"
                "https://alphafold.ebi.ac.uk/files/AF-A0A2M9X344-F1-model_v6.pdb\n"
                "Q1R415\tAFDB:AF-Q1R415-F1\t0.98\t1.0\t1.0\t60\t60\t88\t"
                "https://alphafold.ebi.ac.uk/files/AF-Q1R415-F1-model_v6.pdb\n"
                "A0A7S9SVH6\tAFDB:AF-A0A2M9X3G8-F1\t0.9965277777777778\t1.0\t1.0\t60\t60\t88\t"
                "https://alphafold.ebi.ac.uk/files/AF-A0A2M9X3G8-F1-model_v6.pdb\n"
                "H9L4R6\tAFDB:AF-H9L4R6-F1\t1.0\t1.0\t1.0\t60\t60\t45\t"
                "https://alphafold.ebi.ac.uk/files/AF-H9L4R6-F1-model_v6.pdb\n",
                encoding="utf-8",
            )
            pdb_hits = root / "pdb_hits.tsv"
            pdb_hits.write_text(
                "protein_id\tmatch\tseq_id\tquery_coverage\ttarget_coverage\tqlen\tslen\n"
                "Q1R415\t2cgj_A\t0.983640081799591\t1.0\t1.0\t60\t60\n"
                "A0A7S9SVH6\t1pdl_A\t1.0\t1.0\t1.0\t60\t60\n",
                encoding="utf-8",
            )
            config = {
                "run": {"work_dir": str(root / "run"), "input_fasta": str(fasta)},
                "decision_tree": {
                    "afdb_hits_tsv": str(afdb_hits),
                    "pdb_hits_tsv": str(pdb_hits),
                    "use_afdb_api": False,
                },
                "structure_prediction": {
                    "enabled_steps": ["GETSEQS", "AFDB", "MODELLER", "LOMETS", "DITASSER"],
                },
            }

            run = prepare_decision_tree_run(config)
            self.assertEqual(run_decision_tree(run), 0)

            decisions = {
                line.split("\t")[0]: line.split("\t")
                for line in run.decisions.read_text(encoding="utf-8").splitlines()
            }
            self.assertEqual(decisions["Q1R393"][1], "AFDB")
            self.assertEqual(decisions["A0A080WMA9"][1], "AFDB")
            self.assertEqual(decisions["A0A075B5P5"][1], "MODELLER")
            self.assertEqual(decisions["P03690"][1], "MODELLER")
            self.assertEqual(decisions["Q1R415"][1], "LOMETS")
            self.assertEqual(decisions["A0A7S9SVH6"][1], "LOMETS")
            self.assertEqual(decisions["A0A7S9SVS1"][1], "DITASSER")
            self.assertEqual(decisions["H9L4R6"][1], "DITASSER")

    def test_choose_long_sequence_fallback(self) -> None:
        method, match, score, reason = choose_structure_method(
            "M" * 1200,
            None,
            None,
            {
                "decision_tree": {
                    "long_sequence_length": 900,
                    "long_sequence_method": "DMFOLD",
                },
                "structure_prediction": {"enabled_steps": ["DITASSER", "DMFOLD"]},
            },
        )

        self.assertEqual(
            (method, match, score, reason),
            ("DMFOLD", "", "0", "long_sequence_fallback"),
        )

    def test_choose_exact_afdb(self) -> None:
        result = choose_structure_method(
            "M" * 100,
            _match("P1", "AFDB:AF-P1-F1", seq_id=1.0, qcov=1.0, tcov=1.0, plddt=90),
            None,
            {"structure_prediction": {"enabled_steps": ["AFDB", "MODELLER", "DITASSER"]}},
        )

        self.assertEqual(result, ("AFDB", "AFDB:AF-P1-F1", "1", "afdb_exact_high_confidence"))

    def test_choose_pdb_before_close_afdb(self) -> None:
        result = choose_structure_method(
            "M" * 100,
            _match("P1", "AFDB:AF-P1-F1", seq_id=0.98, qcov=1.0, tcov=1.0, plddt=90),
            _match("P1", "1abc_A", seq_id=1.0, qcov=1.0, tcov=1.0),
            {"structure_prediction": {"enabled_steps": ["AFDB", "MODELLER", "LOMETS", "DITASSER"]}},
        )

        self.assertEqual(result, ("LOMETS", "1abc_A", "1", "pdb_high_identity"))

    def test_choose_close_afdb_template(self) -> None:
        result = choose_structure_method(
            "M" * 100,
            _match("P1", "AFDB:AF-P1-F1", seq_id=0.96, qcov=1.0, tcov=1.0, plddt=90),
            None,
            {"structure_prediction": {"enabled_steps": ["AFDB", "MODELLER", "DITASSER"]}},
        )

        self.assertEqual(
            result,
            ("MODELLER", "AFDB:AF-P1-F1", "0.96", "afdb_template_high_confidence"),
        )

    def test_choose_ditasser_for_low_quality_or_missing_hits(self) -> None:
        low_quality_afdb = choose_structure_method(
            "M" * 100,
            _match("P1", "AFDB:AF-P1-F1", seq_id=1.0, qcov=1.0, tcov=1.0, plddt=45),
            None,
            {"structure_prediction": {"enabled_steps": ["AFDB", "MODELLER", "DITASSER"]}},
        )
        no_hits = choose_structure_method(
            "M" * 100,
            None,
            None,
            {"structure_prediction": {"enabled_steps": ["DITASSER"]}},
        )

        self.assertEqual(low_quality_afdb, ("DITASSER", "", "0", "no_high_confidence_template"))
        self.assertEqual(no_hits, ("DITASSER", "", "0", "no_high_confidence_template"))

    def test_choose_dmfold_when_ditasser_disabled(self) -> None:
        result = choose_structure_method(
            "M" * 100,
            None,
            None,
            {"structure_prediction": {"enabled_steps": ["DMFOLD"]}},
        )

        self.assertEqual(result, ("DMFOLD", "", "0", "fallback"))

    def test_choose_needs_user_model_when_no_method_enabled(self) -> None:
        result = choose_structure_method(
            "M" * 100,
            None,
            None,
            {"structure_prediction": {"enabled_steps": []}},
        )

        self.assertEqual(result, ("NEEDS_USER_MODEL", "", "0", "no_enabled_structure_method"))


class _FakeResponse:
    def __init__(self, payload) -> None:
        import json

        self.payload = json.dumps(payload).encode("utf-8")

    def __enter__(self):
        return self

    def __exit__(self, _exc_type, _exc_value, _traceback) -> None:
        return None

    def read(self) -> bytes:
        return self.payload


def _match(
    protein_id: str,
    match: str,
    seq_id: float,
    qcov: float,
    tcov: float,
    plddt=None,
) -> Match:
    return Match(
        protein_id=protein_id,
        match=match,
        evalue=0.0,
        seq_id=seq_id,
        query_coverage=qcov,
        target_coverage=tcov,
        query_length=100,
        target_length=100,
        align_length=100,
        plddt=plddt,
    )


if __name__ == "__main__":
    unittest.main()
