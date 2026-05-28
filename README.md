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

This repository provides wrapper code, a native decision tree, public
configuration templates, tests, and small examples. It does not redistribute
large databases, model files, containers, or heavyweight third-party
structure/function software. Users are responsible for installing those
external tools and databases, then providing their install paths in a YAML
config file.

Main third-party components used by the workflow include DIAMOND,
PDB/AlphaFold Database resources, supported structure-prediction software,
Singularity or Apptainer, and StarFunc. See `THIRD_PARTY.md` for details.

## Installation

Clone the repository and create the recommended Conda environment:

```bash
git clone git@github.com:QuanEvans/UTI89-Structure-Function-Pipeline.git
cd UTI89-Structure-Function-Pipeline
conda env create -f environment.yml
conda activate uti89-structure-function
pip install -e .
```

The same environment can run the wrapper and the native decision tree. If you
want a smaller environment only for the decision tree and PDB search helper,
use:

```bash
conda env create -f config/decision_tree.environment.yaml
conda activate uti89-decision-tree
```

## Quick Start

Copy the public config template:

```bash
cp config/template.yaml config/my_run.yaml
```

Edit these required paths in `config/my_run.yaml`:

- `run.work_dir`: output directory for this analysis.
- `run.input_fasta`: input protein FASTA.
- `structure_prediction.tools.*`: paths to Singularity/Apptainer structure
  containers and structure-prediction installs.
- `function_prediction.starfunc_sif`: StarFunc container.
- `function_prediction.starfunc_database`: StarFunc database directory.

The bundled decision tree is used by default. It can consume precomputed AFDB
or PDB hit reports, query the AlphaFold DB API by UniProt accession, query the
AlphaFold DB sequence-summary endpoint with the input protein sequence, and run
DIAMOND against a local PDB SEQRES database.

To build the optional PDB SEQRES DIAMOND database:

```bash
python3 scripts/setup_pdb_diamond.py --out-dir /path/to/pdb-diamond
```

Then set:

```yaml
decision_tree:
  pdb_dmnd: /path/to/pdb-diamond/pdb.dmnd
  diamond_executable: /path/to/diamond
```

If `diamond_executable` is omitted, `diamond` must be available on `PATH`.

AFDB hit reports can include a `pdb_path` column with either a local path or an
HTTP(S) AFDB model URL. The pipeline downloads current AFDB model files and
uses Gemmi to convert CIF/mmCIF models to PDB before passing them to downstream
structure steps.

For a local workstation, keep:

```yaml
execution:
  backend: local
  local_cores: 8
  use_modules: false
```

For a Slurm cluster, set:

```yaml
execution:
  backend: slurm
  use_modules: true

slurm:
  account: your_account
  partition: your_partition
```

The `slurm:` config section is only required when `execution.backend: slurm`.
For local runs, omit Slurm settings.

Run the full pipeline and submit generated jobs:

```bash
python3 scripts/run_pipeline.py --config config/my_run.yaml --submit --wait
```

Without `--submit`, the command prepares decision, structure, and function
scripts but does not execute the generated structure/function jobs:

```bash
python3 scripts/run_pipeline.py --config config/my_run.yaml
```

For local execution, `--submit` runs generated scripts with Bash. For Slurm
execution, `--submit` submits them with `sbatch`.

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
