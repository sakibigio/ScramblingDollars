#!/bin/bash
# Markov-estimation margin experiments on the committed (cbase) sigma series:
#   full_rho99    : full sample, stationarity screen rho < 0.99
#   end2021       : estimate through Dec 2021 (row 251) — includes COVID spike
#   end2021_rho99 : same window + stationarity screen
# Live files backed up/restored; outputs tagged data/MS_params_ms_<tag>.csv
set -e
cd "$(dirname "$0")"
JULIA="$HOME/.juliaup/bin/julia"
BK="_ms_margins_baseline"; mkdir -p "$BK/data"
FILES=(RW_shock.csv data/MS_params.csv data/MS_sigma_us_prob.csv data/MS_sigma_us_counterfactuals.csv data/MS_sigma_us_simul.csv data/MS_sigma_eu_prob.csv data/MS_sigma_eu_counterfactuals.csv data/MS_sigma_eu_simul.csv)
for f in "${FILES[@]}"; do [ -f "$f" ] && cp "$f" "$BK/$f"; done
restore() { for f in "${FILES[@]}"; do if [ -f "$BK/$f" ]; then cp "$BK/$f" "$f"; fi; done; echo "[restored live files]"; }
trap restore EXIT
cp RW_shock_cbase.csv RW_shock.csv

run_cfg () {
    TAG=$1; shift
    echo "=== $TAG ($*) ==="
    env "$@" "$JULIA" markov_estimation.jl > "ms_margin_${TAG}_log.txt" 2>&1 || { echo "  FAILED (see log)"; return 0; }
    cp data/MS_params.csv "data/MS_params_ms_${TAG}.csv"
    cp data/MS_sigma_us_prob.csv "data/MS_sigma_us_prob_ms_${TAG}.csv"
    grep -E "mu_sigma_us|rho_sigma_us|^P" "data/MS_params_ms_${TAG}.csv" | head -6
}
run_cfg full_rho99    MS_EST_END=0   MS_RHO_MAX=0.99
run_cfg end2021       MS_EST_END=251
run_cfg end2021_rho99 MS_EST_END=251 MS_RHO_MAX=0.99
echo "=== DONE ==="
