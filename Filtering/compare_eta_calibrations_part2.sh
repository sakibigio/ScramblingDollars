#!/bin/bash
# Part 2 of eta sweep: covers eta = 0.35 and 0.25 only.
# Uses nohup-safe invocation so OS session changes don't kill it.
# IMPORTANT: does NOT abort on a single matlab failure — runs all etas.

cd "$(dirname "$0")"

ETA_VALUES=(0.35 0.25)
MATLAB="${MATLAB_BIN:-/Applications/MATLAB_R2025b.app/bin/matlab}"
JULIA="$(which julia)"
OVERLEAF="${SCRAMBLING_QUANTFIGS:-/Users/sakibigio/Dropbox/Apps/Overleaf/ScramblingDollarsLiquidity_NewVersion_Restud/quantfigs}"

for ETA in "${ETA_VALUES[@]}"; do
    TAG="eta$(echo "$ETA" | awk '{printf "%.0f", $1*100}')"
    echo
    echo "==================================================================="
    echo "  η = $ETA  (tag: $TAG)  -- $(date +%H:%M:%S)"
    echo "==================================================================="

    # Write override
    "$MATLAB" -batch "eta_override_val=$ETA; save('_eta_session_override.mat','eta_override_val');"
    if [ $? -ne 0 ]; then echo "  override write FAILED"; continue; fi

    # Filter
    echo "[1/3] main_filter.m  $(date +%H:%M:%S)"
    "$MATLAB" -batch "main_filter; quit force"
    if [ $? -ne 0 ]; then echo "  main_filter FAILED for $TAG"; continue; fi
    cp RW_shock.csv "RW_shock_${TAG}.csv"

    # Julia Markov
    echo "[2/3] markov_estimation.jl  $(date +%H:%M:%S)"
    "$JULIA" markov_estimation.jl
    if [ $? -ne 0 ]; then echo "  julia FAILED for $TAG"; continue; fi
    cp data/MS_params.csv         "data/MS_params_${TAG}.csv"
    cp data/MS_sigma_us_prob.csv  "data/MS_sigma_us_prob_${TAG}.csv"

    # main_LFX
    echo "[3/3] main_LFX.m  $(date +%H:%M:%S)"
    "$MATLAB" -batch "matching_type=1; save('temp_matching_type.mat','matching_type'); main_LFX; quit force"
    if [ $? -ne 0 ]; then echo "  main_LFX FAILED for $TAG"; continue; fi
    cp "$OVERLEAF/Mod_Moments_cd.tex"           "$OVERLEAF/Mod_Moments_cd_${TAG}.tex"
    cp "$OVERLEAF/Mod_CIP_Moments_cd.tex"       "$OVERLEAF/Mod_CIP_Moments_cd_${TAG}.tex"
    cp "$OVERLEAF/Mod_SwitchDev_Moments_cd.tex" "$OVERLEAF/Mod_SwitchDev_Moments_cd_${TAG}.tex"
    echo "  ✓ saved $TAG outputs  $(date +%H:%M:%S)"
done

rm -f _eta_session_override.mat
echo
echo "DONE $(date +%H:%M:%S)"
