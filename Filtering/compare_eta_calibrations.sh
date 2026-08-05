#!/bin/bash
#
# compare_eta_calibrations.sh
#
# Sweep the full pipeline (main_filter → markov_estimation.jl → main_LFX) over
# several values of eta and tag the outputs for moment comparison.
#
# Hooks:
#   - params.m honors `_eta_session_override.mat` (written by this script) to
#     temporarily override η for both the filter and the model.
#   - markov_estimation.jl is unchanged — its inputs (filtered σ series) reflect
#     the chosen η automatically.
#
# Outputs per eta value (suffix _etaXX where XX = eta * 100, e.g. eta50):
#   data/MS_params_etaXX.csv             (Markov parameters)
#   data/MS_sigma_us_prob_etaXX.csv      (Hamilton-filtered probabilities)
#   RW_shock_etaXX.csv                   (filtered σ series)
#   <overleaf>/Mod_*_cd_etaXX.tex        (moment tables)
#
# Usage: ./compare_eta_calibrations.sh
#

set -e
cd "$(dirname "$0")"

# ---- Configuration ---------------------------------------------------------
ETA_VALUES=(0.50 0.35 0.25)
MATLAB="${MATLAB_BIN:-/Applications/MATLAB_R2025b.app/bin/matlab}"
JULIA="$(which julia)"
OVERLEAF="${SCRAMBLING_QUANTFIGS:-/Users/sakibigio/Dropbox/Apps/Overleaf/ScramblingDollarsLiquidity_NewVersion_Restud/quantfigs}"

# Sanity checks
if [ ! -x "$MATLAB" ]; then
    echo "ERROR: MATLAB not found at $MATLAB"
    echo "Set MATLAB_BIN=/path/to/matlab, or edit the MATLAB= line above."
    exit 1
fi
if [ -z "$JULIA" ]; then
    echo "ERROR: julia not on PATH."
    exit 1
fi

# ---- Backup baseline state -------------------------------------------------
BACKUP_DIR="_eta_sensitivity_baseline"
mkdir -p "$BACKUP_DIR"
[ -f data/MS_params.csv ]            && cp data/MS_params.csv            "$BACKUP_DIR/" || true
[ -f data/MS_sigma_us_prob.csv ]     && cp data/MS_sigma_us_prob.csv     "$BACKUP_DIR/" || true
[ -f RW_shock.csv ]                  && cp RW_shock.csv                  "$BACKUP_DIR/" || true
[ -f "$OVERLEAF/Mod_Moments_cd.tex" ]      && cp "$OVERLEAF/Mod_Moments_cd.tex"      "$BACKUP_DIR/" || true
[ -f "$OVERLEAF/Mod_CIP_Moments_cd.tex" ]  && cp "$OVERLEAF/Mod_CIP_Moments_cd.tex"  "$BACKUP_DIR/" || true
[ -f "$OVERLEAF/Mod_SwitchDev_Moments_cd.tex" ] && cp "$OVERLEAF/Mod_SwitchDev_Moments_cd.tex" "$BACKUP_DIR/" || true
echo "Baseline state backed up to $BACKUP_DIR/"
echo

# ---- Eta sweep -------------------------------------------------------------
for ETA in "${ETA_VALUES[@]}"; do
    # Tag: eta50, eta35, eta25
    TAG="eta$(echo "$ETA" | awk '{printf "%.0f", $1*100}')"
    echo "==================================================================="
    echo "  ETA SWEEP — η = $ETA  (tag: $TAG)"
    echo "==================================================================="

    # 1) Write override file
    "$MATLAB" -batch "eta_override_val=$ETA; save('_eta_session_override.mat','eta_override_val'); fprintf('  override set: eta=%.3f\n', eta_override_val);"

    # 2) Run filter (writes RW_shock.csv)
    echo "[1/3] main_filter.m..."
    "$MATLAB" -batch "main_filter; quit force"
    cp RW_shock.csv "RW_shock_${TAG}.csv"

    # 3) Run Markov estimation (reads RW_shock.csv, writes data/MS_params.csv)
    echo "[2/3] markov_estimation.jl..."
    "$JULIA" markov_estimation.jl
    cp data/MS_params.csv         "data/MS_params_${TAG}.csv"
    cp data/MS_sigma_us_prob.csv  "data/MS_sigma_us_prob_${TAG}.csv"

    # 4) Run main_LFX (reads MS_params.csv, writes Mod_*_cd.tex on Overleaf)
    echo "[3/3] main_LFX.m..."
    "$MATLAB" -batch "matching_type=1; save('temp_matching_type.mat','matching_type'); main_LFX; quit force"
    cp "$OVERLEAF/Mod_Moments_cd.tex"           "$OVERLEAF/Mod_Moments_cd_${TAG}.tex"
    cp "$OVERLEAF/Mod_CIP_Moments_cd.tex"       "$OVERLEAF/Mod_CIP_Moments_cd_${TAG}.tex"
    cp "$OVERLEAF/Mod_SwitchDev_Moments_cd.tex" "$OVERLEAF/Mod_SwitchDev_Moments_cd_${TAG}.tex"

    echo "  -> saved tagged outputs for $TAG"
    echo
done

# ---- Cleanup ---------------------------------------------------------------
rm -f _eta_session_override.mat

# Restore baseline live state so subsequent ad-hoc runs use the original config
[ -f "$BACKUP_DIR/MS_params.csv" ]           && cp "$BACKUP_DIR/MS_params.csv"           data/MS_params.csv
[ -f "$BACKUP_DIR/MS_sigma_us_prob.csv" ]    && cp "$BACKUP_DIR/MS_sigma_us_prob.csv"    data/MS_sigma_us_prob.csv
[ -f "$BACKUP_DIR/RW_shock.csv" ]            && cp "$BACKUP_DIR/RW_shock.csv"            RW_shock.csv
echo "Live state restored from $BACKUP_DIR/"

# ---- Summary ---------------------------------------------------------------
echo
echo "==================================================================="
echo "  DONE.  Tagged outputs:"
for ETA in "${ETA_VALUES[@]}"; do
    TAG="eta$(echo "$ETA" | awk '{printf "%.0f", $1*100}')"
    echo "    η = $ETA  ($TAG)"
    echo "      data/MS_params_${TAG}.csv"
    echo "      RW_shock_${TAG}.csv"
    echo "      $OVERLEAF/Mod_Moments_cd_${TAG}.tex"
    echo "      $OVERLEAF/Mod_CIP_Moments_cd_${TAG}.tex"
done
echo "==================================================================="
