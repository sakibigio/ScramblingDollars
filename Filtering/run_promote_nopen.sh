#!/bin/bash
# PROMOTION of the nopen_lcr baseline (lambda=1.392, iota=0.0951, eta=0.632).
# Live pipeline + Overleaf tables move to the new calibration; dated backups kept.
set -e
cd "$(dirname "$0")"
MATLAB="${MATLAB_BIN:-/Applications/MATLAB_R2025b.app/bin/matlab}"
JULIA="$HOME/.juliaup/bin/julia"
OVERLEAF="${SCRAMBLING_QUANTFIGS:-/Users/sakibigio/Dropbox/Apps/Overleaf/ScramblingDollarsLiquidity_NewVersion_Restud/quantfigs}"
STAMP=prepromote_20260731

echo "=== [0] Backups ==="
cp _calibration_override.mat "_calibration_override_${STAMP}.mat"
cp RW_shock.csv "RW_shock_${STAMP}.csv"
cp data/MS_params.csv "data/MS_params_${STAMP}.csv"
for f in Mod_Moments_cd Mod_CIP_Moments_cd Mod_SwitchDev_Moments_cd Mod_Moments Mod_CIP_Moments Mod_SwitchDev_Moments Data_CIP_Moments Data_Moments tab_markov_estimates_cd; do
    [ -f "$OVERLEAF/$f.tex" ] && cp "$OVERLEAF/$f.tex" "$OVERLEAF/${f}_${STAMP}.tex"
done
echo "backups stamped _${STAMP}"

echo "=== [1] Promote calibration + warm start ==="
cp _calibration_override_nopen_lcr.mat _calibration_override.mat
cp data/initguess_nopen_lcr.mat data/initguess.mat

echo "=== [2] main_filter (promoted calibration; LCR mu; writes RW_shock + filter figures) ==="
"$MATLAB" -batch "main_filter; quit force" > promote_filter_log.txt 2>&1
python3 - <<'PY'
import csv
a=[float(r['sigma_us']) for r in csv.DictReader(open('RW_shock.csv'))]
b=[float(r['sigma_us']) for r in csv.DictReader(open('RW_shock_est_nopen_lcr.csv'))]
d=max(abs(x-y) for x,y in zip(a,b))
print('RW_shock vs tagged nopen series: max |diff| = %.2e over %d rows' % (d, len(a)))
PY

echo "=== [3] Markov estimation (full MS_* set) ==="
"$JULIA" markov_estimation.jl > promote_markov_log.txt 2>&1
grep -E "mu_sigma_us_r|P11|P22" data/MS_params.csv

echo "=== [4] Markov table + data-moment tables ==="
"$MATLAB" -batch "regen_markov_tex; quit force" > promote_tex_log.txt 2>&1 || echo "  regen_markov_tex FAILED (check log)"
"$MATLAB" -batch "compute_data_moments; quit force" > promote_datamom_log.txt 2>&1 || echo "  compute_data_moments FAILED (check log)"

echo "=== [5] Mirror model moment tables (from the nopen_lcr main_LFX pass) ==="
for t in Mod_Moments_cd Mod_CIP_Moments_cd Mod_SwitchDev_Moments_cd; do
    cp "$OVERLEAF/${t}_nopen_lcr.tex" "$OVERLEAF/${t}.tex"
    base=${t%_cd}
    cp "$OVERLEAF/${t}_nopen_lcr.tex" "$OVERLEAF/${base}.tex"
    echo "  -> $t.tex and $base.tex"
done

echo "=== [6] Regenerated data tables ==="
cat "$OVERLEAF/Data_CIP_Moments.tex" 2>/dev/null
echo "=== PROMOTION DONE ==="
