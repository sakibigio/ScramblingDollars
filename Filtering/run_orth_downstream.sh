#!/bin/bash
#
# run_orth_downstream.sh
#
# Downstream Markov check for the alternative identification-penalty
# calibrations (orth_mode = trend, kpss). For each mode:
#   1) swap _calibration_override_<mode>.mat into _calibration_override.mat
#      (params.m honors it for both filter and model),
#   2) run main_filter.m  -> RW_shock.csv          (tagged RW_shock_<mode>.csv)
#   3) run markov_estimation.jl -> data/MS_params.csv etc. (tagged _<mode>)
# All live pipeline files are backed up first and restored at the end
# (trap ensures restore even on failure), so the committed baseline state
# is untouched.
#
# Usage: ./run_orth_downstream.sh

set -e
cd "$(dirname "$0")"

MATLAB="${MATLAB_BIN:-/Applications/MATLAB_R2025b.app/bin/matlab}"
JULIA="$HOME/.juliaup/bin/julia"
MODES=(${ORTH_MODES:-trend kpss})   # override with e.g. ORTH_MODES=eta50

BACKUP_DIR="_orth_downstream_baseline"
mkdir -p "$BACKUP_DIR"

LIVE_FILES=(
    "_calibration_override.mat"
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

for MODE in "${MODES[@]}"; do
    echo
    echo "==================================================================="
    echo "  ORTH DOWNSTREAM — mode = $MODE"
    echo "==================================================================="
    if [ ! -f "_calibration_override_${MODE}.mat" ]; then
        echo "ERROR: _calibration_override_${MODE}.mat not found; skipping."
        continue
    fi
    cp "_calibration_override_${MODE}.mat" "_calibration_override.mat"

    echo "[1/2] main_filter.m ..."
    "$MATLAB" -batch "main_filter; quit force" > "orth_${MODE}_filter_log.txt" 2>&1
    cp RW_shock.csv "RW_shock_${MODE}.csv"
    echo "      -> RW_shock_${MODE}.csv"

    echo "[2/2] markov_estimation.jl ..."
    "$JULIA" markov_estimation.jl > "orth_${MODE}_markov_log.txt" 2>&1
    cp data/MS_params.csv        "data/MS_params_${MODE}.csv"
    cp data/MS_sigma_us_prob.csv "data/MS_sigma_us_prob_${MODE}.csv"
    [ -f data/MS_params_eu.csv ] && cp data/MS_params_eu.csv "data/MS_params_eu_${MODE}.csv"
    echo "      -> data/MS_params_${MODE}.csv"
    echo
    echo "--- MS_params ($MODE) ---"
    cat "data/MS_params_${MODE}.csv"
done

echo
echo "=== DONE ==="
