# Third-Party Software and Data

This repository provides wrappers and configuration templates. It does not
redistribute the external software, containers, or databases required by the
pipeline.

Users are responsible for installing, licensing, and citing the tools and data
they configure, including:

- Decision-tree source workflow and Snakemake environment.
- DIAMOND databases and FASTA sequence libraries.
- AlphaFold DB foldcomp data.
- Structure tools such as D-I-TASSER, LOMETS, MODELLER, and optional DMFold.
- Singularity or Apptainer.
- StarFunc container and StarFunc database.

Copy `config/template.yaml` or `config/user_models.example.yaml` and replace
all paths with your own installs.
