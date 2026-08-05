#!/bin/bash
# 2x3 matrix: {raw mu, LCR mu} x {corr, trend, kpss}. Raw-mu estimations reuse
# existing results; LCR-mu estimations run fresh. Then coherent sigma series
# per cell + truncated Markov. Live files backed up/restored.
set -e
cd "$(dirname "$0")"
MATLAB="${MATLAB_BIN:-/Applications/MATLAB_R2025b.app/bin/matlab}"
JULIA="$HOME/.juliaup/bin/julia"

echo "=== [1/3] LCR-mu estimations (corr, trend, kpss) ==="
for MODE in corr trend kpss; do
    echo "--- estimating: $MODE (LCR mu) ---"
    ORTH_MODE=$MODE MU_LCR=1 "$MATLAB" -batch run_estimate_orth > "est_${MODE}_lcr_log.txt" 2>&1 \
        || { echo "  FAILED (see est_${MODE}_lcr_log.txt)"; continue; }
    grep -E "lambda    =|iota_ss   =|eta       =" "est_${MODE}_lcr_log.txt" | head -3
done

echo "=== [2/3] Coherent sigma series for all 6 cells ==="
"$MATLAB" -batch gen_est_series > gen_est_series_log.txt 2>&1
grep -E "RW_shock_est|SKIP" gen_est_series_log.txt

echo "=== [3/3] Truncated Markov per cell ==="
BK="_muscore_baseline"; mkdir -p "$BK/data"
FILES=(RW_shock.csv data/MS_params.csv data/MS_sigma_us_prob.csv data/MS_sigma_us_counterfactuals.csv data/MS_sigma_us_simul.csv data/MS_sigma_eu_prob.csv data/MS_sigma_eu_counterfactuals.csv data/MS_sigma_eu_simul.csv)
for f in "${FILES[@]}"; do [ -f "$f" ] && cp "$f" "$BK/$f"; done
restore() { for f in "${FILES[@]}"; do if [ -f "$BK/$f" ]; then cp "$BK/$f" "$f"; fi; done; echo "[restored]"; }
trap restore EXIT
for CELL in est_corr est_trend est_kpss est_corr_lcr est_trend_lcr est_kpss_lcr; do
    [ -f "RW_shock_${CELL}.csv" ] || { echo "$CELL: no series, skipping"; continue; }
    cp "RW_shock_${CELL}.csv" RW_shock.csv
    "$JULIA" markov_estimation.jl > "markov_${CELL}_log.txt" 2>&1 || { echo "$CELL: MARKOV FAILED"; continue; }
    cp data/MS_params.csv "data/MS_params_${CELL}.csv"
    echo "--- $CELL ---"
    grep -E "mu_sigma_us|Sigma_sigma_us|^P11|^P22" "data/MS_params_${CELL}.csv"
done
echo "=== DONE ==="
