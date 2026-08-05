#!/bin/bash
#
# run_dual_lfx.sh
#
# Clean dual pipeline pass: for each candidate calibration (cbase = committed
# lambda=1.835/eta=0.611, eta50 = fixed eta=0.5, lambda=1.248/iota=0.070),
# stage its tagged filter output, regenerate ALL Markov files consistently,
# run main_LFX, and tag the Overleaf moment tables. Live pipeline files and
# the untagged Overleaf tables are backed up first and restored at the end,
# so nothing untagged changes until a winner is promoted.
#
# PRECONDITIONS:
#   - No robustness sweeps or other MATLAB/Julia writers running.
#   - The nonlinear-model warm-start work (other session) is merged; if it
#     introduced cache files for main_LFX, list them in CACHE_FILES below so
#     each tag's run is staged/restored cleanly.
#
# Usage: ./run_dual_lfx.sh

set -e
cd "$(dirname "$0")"

MATLAB="${MATLAB_BIN:-/Applications/MATLAB_R2025b.app/bin/matlab}"
JULIA="$HOME/.juliaup/bin/julia"
. ./overleaf_dir.sh            # sets $OVERLEAF (repo-local unless this machine has the sync folder)
TAGS=(${LFX_TAGS:-cbase eta50})

# Warm-start cache: solve_global.m loads data/initguess.mat (with a size-
# mismatch guard) and re-saves it after every successful solve. Tagged per
# calibration below so promotion can restore the winner's solution.
CACHE_FILES=("data/initguess.mat")

BACKUP_DIR="_dual_lfx_baseline"
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
    "data/initguess.mat"
)
TEX_FILES=(
    "Mod_Moments_cd.tex"
    "Mod_CIP_Moments_cd.tex"
    "Mod_SwitchDev_Moments_cd.tex"
)

echo "=== Backing up live pipeline + Overleaf tables to $BACKUP_DIR/ ==="
for f in "${LIVE_FILES[@]}"; do
    if [ -f "$f" ]; then
        mkdir -p "$BACKUP_DIR/$(dirname "$f")"
        cp "$f" "$BACKUP_DIR/$f"; echo "  backed up: $f"
    fi
done
mkdir -p "$BACKUP_DIR/overleaf"
for f in "${TEX_FILES[@]}"; do
    [ -f "$OVERLEAF/$f" ] && { cp "$OVERLEAF/$f" "$BACKUP_DIR/overleaf/$f"; echo "  backed up: overleaf/$f"; }
done

restore_live() {
    echo "=== Restoring live pipeline + Overleaf tables ==="
    for f in "${LIVE_FILES[@]}"; do
        if [ -f "$BACKUP_DIR/$f" ]; then cp "$BACKUP_DIR/$f" "$f"; echo "  restored: $f"; fi
    done
    for f in "${TEX_FILES[@]}"; do
        if [ -f "$BACKUP_DIR/overleaf/$f" ]; then cp "$BACKUP_DIR/overleaf/$f" "$OVERLEAF/$f"; echo "  restored: overleaf/$f"; fi
    done
}
trap restore_live EXIT

for TAG in "${TAGS[@]}"; do
    echo
    echo "==================================================================="
    echo "  DUAL LFX PASS — calibration = $TAG"
    echo "==================================================================="
    if [ ! -f "_calibration_override_${TAG}.mat" ] || [ ! -f "RW_shock_${TAG}.csv" ]; then
        echo "ERROR: tagged inputs for $TAG missing; skipping."
        continue
    fi

    cp "_calibration_override_${TAG}.mat" "_calibration_override.mat"
    cp "RW_shock_${TAG}.csv" "RW_shock.csv"

    echo "[1/2] markov_estimation.jl (regenerate ALL MS_* consistently) ..."
    "$JULIA" markov_estimation.jl > "dual_${TAG}_markov_log.txt" 2>&1

    echo "[2/2] main_LFX.m ..."
    "$MATLAB" -batch "matching_type=1; save('temp_matching_type.mat','matching_type'); main_LFX; quit force" \
        > "dual_${TAG}_lfx_log.txt" 2>&1

    for f in "${TEX_FILES[@]}"; do
        [ -f "$OVERLEAF/$f" ] && { cp "$OVERLEAF/$f" "$OVERLEAF/${f%.tex}_${TAG}.tex"; echo "  -> ${f%.tex}_${TAG}.tex"; }
    done
    for f in "${CACHE_FILES[@]}"; do
        [ -f "$f" ] && cp "$f" "${f%.mat}_${TAG}.mat"
    done
done

echo
echo "=== DONE — compare tagged tables, then promote the winner ==="
