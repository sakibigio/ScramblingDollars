#!/bin/bash
#
# run_orth_fullsample.sh
#
# Full-sample (no pre-2020 truncation) Markov estimation for the baseline
# (committed corr calibration) and the trend/kpss alternative calibrations.
# Reuses the already-filtered sigma series:
#   baseline: _orth_downstream_baseline/RW_shock.csv (committed pipeline output)
#   trend:    RW_shock_trend.csv   (from run_orth_downstream.sh)
#   kpss:     RW_shock_kpss.csv    (from run_orth_downstream.sh)
# Runs markov_estimation.jl with MS_EST_END=0 (full sample) and tags outputs
# MS_params_<tag>_full.csv / MS_sigma_us_prob_<tag>_full.csv.
# Live pipeline files are backed up and restored via trap.
#
# Usage: ./run_orth_fullsample.sh

set -e
cd "$(dirname "$0")"

JULIA="$HOME/.juliaup/bin/julia"

BACKUP_DIR="_orth_fullsample_baseline"
mkdir -p "$BACKUP_DIR"

LIVE_FILES=(
    "RW_shock.csv"
    "data/MS_params.csv"
    "data/MS_sigma_us_prob.csv"
    "data/MS_sigma_us_counterfactuals.csv"
    "data/MS_sigma_us_simul.csv"
    "data/MS_sigma_eu_prob.csv"
    "data/MS_sigma_eu_counterfactuals.csv"
    "data/MS_sigma_eu_simul.csv"
    "data/MS_params_eu.csv"
)

echo "=== Backing up live pipeline files to $BACKUP_DIR/ ==="
for f in "${LIVE_FILES[@]}"; do
    if [ -f "$f" ]; then
        mkdir -p "$BACKUP_DIR/$(dirname "$f")"
        cp "$f" "$BACKUP_DIR/$f"
        echo "  backed up: $f"
    fi
done

restore_live() {
    echo "=== Restoring live pipeline files from $BACKUP_DIR/ ==="
    for f in "${LIVE_FILES[@]}"; do
        if [ -f "$BACKUP_DIR/$f" ]; then
            cp "$BACKUP_DIR/$f" "$f"
            echo "  restored: $f"
        fi
    done
}
trap restore_live EXIT

for TAG in ${ORTH_TAGS:-baseline trend kpss}; do
    case "$TAG" in
        baseline) SRC="_orth_downstream_baseline/RW_shock.csv" ;;
        *)        SRC="RW_shock_${TAG}.csv" ;;
    esac
    echo
    echo "==================================================================="
    echo "  FULL-SAMPLE MARKOV — calibration = $TAG  (shock: $SRC)"
    echo "==================================================================="
    if [ ! -f "$SRC" ]; then
        echo "ERROR: $SRC not found; skipping."
        continue
    fi
    cp "$SRC" RW_shock.csv

    MS_EST_END=0 "$JULIA" markov_estimation.jl > "orth_${TAG}_markov_full_log.txt" 2>&1
    cp data/MS_params.csv        "data/MS_params_${TAG}_full.csv"
    cp data/MS_sigma_us_prob.csv "data/MS_sigma_us_prob_${TAG}_full.csv"
    echo "      -> data/MS_params_${TAG}_full.csv"
    echo
    echo "--- MS_params ($TAG, FULL sample) ---"
    cat "data/MS_params_${TAG}_full.csv"
done

echo
echo "=== DONE ==="
