# UTI89 Structure-Function Pipeline

This repository contains the workflow used to generate proteome-scale
structure and function predictions for uropathogenic *Escherichia coli*
UTI89. The pipeline starts from protein FASTA sequences, selects or predicts a
high-confidence structural model, and runs StarFunc to assign protein function
annotations.

The accompanying UTI89 structure-function database is available here:

https://seq2fun.dcmb.med.umich.edu/UTI89

## What This Pipeline Does

- Runs a decision tree that chooses an appropriate structure source or modeling
  method for each input protein.
- Uses existing high-confidence PDB or AlphaFold Database models when
  available.
- Prepares structure-prediction jobs for supported modeling methods when no
  suitable existing model is available.
- Accepts user-provided PDB models for users who already have structures.
- Runs StarFunc on staged `model1.pdb` structures to produce function
  prediction outputs.
- Supports direct local Bash execution or Slurm submission through YAML
  configuration.

## Repository Scope

This repository provides wrapper code, public configuration templates, tests,
and small examples. It does not redistribute large databases, model files,
containers, or third-party software. Users are responsible for installing the
required tools and databases, then providing their install paths in a YAML
config file.

Main third-party components used by the workflow include Snakemake, DIAMOND,
Foldcomp, PDB/AlphaFold Database resources, supported structure-prediction
software, Singularity or Apptainer, and StarFunc. See `THIRD_PARTY.md` for
details.

## Quick Start

Copy and edit the public config template:

```bash
cp config/template.yaml config/my_run.yaml
```

Then run the full pipeline:

```bash
python3 scripts/run_pipeline.py --config config/my_run.yaml --submit --wait
```

For a local workstation, set:

```yaml
execution:
  backend: local
  local_cores: 8
  use_modules: false
```

For a Slurm cluster, set `execution.backend: slurm` and add the required
`slurm:` account and partition settings in the config.

## Using Existing Structures

If you already have protein structures, start from:

```bash
cp config/user_models.example.yaml config/my_models.yaml
```

Configure `structure_prediction.user_models` with either a model directory and
filename pattern or a TSV mapping file. In this mode the pipeline skips the
decision tree, stages each structure as
`structures/<protein_id>/model1.pdb`, and runs StarFunc.

## Common Commands

Run only the decision-tree stage:

```bash
python3 scripts/run_pipeline.py --config config/my_run.yaml --stop-after decision
```

Run only structure preparation and submission:

```bash
python3 scripts/run_pipeline.py --config config/my_run.yaml --start-at structure --stop-after structure --submit --wait
```

Run only StarFunc function prediction from staged structures:

```bash
python3 scripts/run_pipeline.py --config config/my_run.yaml --start-at function --submit --wait
```

Run non-submitting smoke checks:

```bash
scripts/smoke_test.sh
```

## Outputs

For a run directory `$RUN_DIR`, expected outputs include:

- `$RUN_DIR/decision_tree_intermediates/decide/decisions.txt`
- `$RUN_DIR/structures/<protein_id>/seq.fasta`
- `$RUN_DIR/structures/<protein_id>/model1.pdb`
- `$RUN_DIR/functions/<protein_id>/consensus.tsv`
- StarFunc component outputs such as `structure.tsv`, `sequence.tsv`,
  `pfam.tsv`, `ppi.tsv`, `naive.tsv`, and `interlabelgo.tsv`

Decision-tree outputs are skipped for own-model runs.

## Citation

If you use this workflow or the UTI89 structure-function database, cite the
associated manuscript and this repository. A BibTeX entry is provided in
`CITATION.bib`.

## License

The wrapper code in this repository is released under the MIT License. Third
party tools and databases are governed by their own licenses and access terms.
