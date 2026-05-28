#!/usr/bin/env bash
set -euo pipefail

CONFIG="${1:-}"
SMOKE_RUN_DIR="${SMOKE_RUN_DIR:-/tmp/uti89_pipeline_smoke_${USER}_$$}"

PYTHONPATH=src python3 -m unittest discover tests
if [ -n "$CONFIG" ]; then
    python3 scripts/run_pipeline.py \
        --config "$CONFIG" \
        --work-dir "$SMOKE_RUN_DIR" \
        --decision-dry-run \
        --stop-after decision
fi

OWN_MODEL_DIR="$SMOKE_RUN_DIR/user_models"
OWN_MODEL_CONFIG="$SMOKE_RUN_DIR/user_models.yaml"
OWN_MODEL_RUN_DIR="$SMOKE_RUN_DIR/user_model_run"
DECISION_CONFIG="$SMOKE_RUN_DIR/decision_tree.yaml"
DECISION_RUN_DIR="$SMOKE_RUN_DIR/decision_tree_run"
mkdir -p "$OWN_MODEL_DIR"
cat > "$DECISION_CONFIG" <<EOF
run:
  name: decision_tree_smoke
  work_dir: $DECISION_RUN_DIR
  input_fasta: $PWD/examples/smoke.fasta
execution:
  backend: local
decision_tree:
  use_afdb_api: false
structure_prediction:
  enabled_steps:
    - DITASSER
EOF
python3 scripts/run_pipeline.py \
    --config "$DECISION_CONFIG" \
    --work-dir "$DECISION_RUN_DIR" \
    --stop-after decision
test -s "$DECISION_RUN_DIR/decision_tree_intermediates/decide/decisions.txt"

printf 'HEADER SMOKE USER MODEL\nEND\n' > "$OWN_MODEL_DIR/myprot.pdb"
cat > "$OWN_MODEL_CONFIG" <<EOF
run:
  name: user_model_smoke
  work_dir: $OWN_MODEL_RUN_DIR
  input_fasta: $PWD/examples/smoke.fasta
execution:
  backend: local
structure_prediction:
  user_models:
    model_dir: $OWN_MODEL_DIR
    pattern: "{protein_id}.pdb"
    link_models: false
function_prediction:
  starfunc_sif: /example/StarFunc.sif
  starfunc_database: /example/database
EOF

python3 scripts/run_pipeline.py \
    --config "$OWN_MODEL_CONFIG" \
    --work-dir "$OWN_MODEL_RUN_DIR" \
    --stop-after function
test -s "$OWN_MODEL_RUN_DIR/structures/myprot/model1.pdb"
test -f "$OWN_MODEL_RUN_DIR/functions/myprot/run_starfunc.sh"

echo "Smoke checks completed using $SMOKE_RUN_DIR"
