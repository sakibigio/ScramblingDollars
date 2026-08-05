#!/bin/bash
# No-penalty, free-eta, LCR-mu estimation (pure moments) + coherent Markov.
set -e
cd "$(dirname "$0")"
MATLAB="${MATLAB_BIN:-/Applications/MATLAB_R2025b.app/bin/matlab}"
JULIA="$HOME/.juliaup/bin/julia"
# w=0 with ORTH_MODE=corr writes _calibration_override_lcr.mat and
# _estimation_result_corr_lcr.mat -- guard the corr_lcr files, then rename.
cp _calibration_override_lcr.mat _cal_corr_lcr_GUARD.mat 2>/dev/null || true
cp _estimation_result_corr_lcr.mat _res_corr_lcr_GUARD.mat 2>/dev/null || true
ORTH_MODE=corr ORTH_W=0 MU_LCR=1 "$MATLAB" -batch run_estimate_orth > est_nopen_lcr_log.txt 2>&1
mv _calibration_override_lcr.mat _calibration_override_nopen_lcr.mat
mv _estimation_result_corr_lcr.mat _estimation_result_nopen_lcr.mat
[ -f _cal_corr_lcr_GUARD.mat ] && mv _cal_corr_lcr_GUARD.mat _calibration_override_lcr.mat
[ -f _res_corr_lcr_GUARD.mat ] && mv _res_corr_lcr_GUARD.mat _estimation_result_corr_lcr.mat
grep -E "lambda    =|iota_ss   =|eta       =|CIP        :|mean\(FF\)   :|mean\(DWS\)  :" est_nopen_lcr_log.txt | head -6

# Coherent series (reuse gen script pattern via result file rename trick)
cp _estimation_result_nopen_lcr.mat _estimation_result_eta50_lcr_TMPSWAP.mat
mv _estimation_result_eta50_lcr.mat _estimation_result_eta50_lcr_REAL.mat
mv _estimation_result_eta50_lcr_TMPSWAP.mat _estimation_result_eta50_lcr.mat
"$MATLAB" -batch gen_eta50_lcr_series > gen_nopen_lcr_log.txt 2>&1
mv RW_shock_est_eta50_lcr.csv RW_shock_est_nopen_lcr.csv
mv _estimation_result_eta50_lcr_REAL.mat _estimation_result_eta50_lcr.mat
grep "RW_shock" gen_nopen_lcr_log.txt

# Markov (truncated) with backup/restore
BK=_nopen_mk; mkdir -p $BK/data
FILES=(RW_shock.csv data/MS_params.csv data/MS_sigma_us_prob.csv data/MS_sigma_us_counterfactuals.csv data/MS_sigma_us_simul.csv data/MS_sigma_eu_prob.csv data/MS_sigma_eu_counterfactuals.csv data/MS_sigma_eu_simul.csv)
for f in "${FILES[@]}"; do [ -f "$f" ] && cp "$f" "$BK/$f"; done
restore() { for f in "${FILES[@]}"; do if [ -f "$BK/$f" ]; then cp "$BK/$f" "$f"; fi; done; echo "[restored]"; }
trap restore EXIT
cp RW_shock_est_nopen_lcr.csv RW_shock.csv
"$JULIA" markov_estimation.jl > markov_nopen_lcr_log.txt 2>&1
cp data/MS_params.csv data/MS_params_est_nopen_lcr.csv
echo "--- nopen LCR Markov (truncated) ---"
grep -E "mu_sigma_us|Sigma_sigma_us|P11|P22" data/MS_params_est_nopen_lcr.csv
