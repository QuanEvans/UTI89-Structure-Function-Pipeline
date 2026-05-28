# Third-Party Software and Data

This repository provides wrappers, configuration templates, and the small
native decision tree needed to select structure-generation routes. It does not
redistribute large external software, containers, or databases required by the
pipeline.

Users are responsible for installing, licensing, and citing the tools and data
they configure, including:

- DIAMOND and a local PDB SEQRES database if PDB template search is enabled.
- AlphaFold DB PDB model files, AlphaFold DB API access for accession or
  sequence lookup, or user-supplied AlphaFold hit reports.
- Structure tools such as D-I-TASSER, LOMETS, MODELLER, and optional DMFold.
- Singularity or Apptainer.
- StarFunc container and StarFunc database.

Copy `config/template.yaml` or `config/user_models.example.yaml` and replace
all paths with your own installs.
