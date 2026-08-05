#!/bin/bash
# Ex-post moments for: committed calibration + Dec-2021 Markov window.
set -e
cd "$(dirname "$0")"
MATLAB="${MATLAB_BIN:-/Applications/MATLAB_R2025b.app/bin/matlab}"
. ./overleaf_dir.sh            # sets $OVERLEAF (repo-local unless this machine has the sync folder)
BK="_end2021_baseline"; mkdir -p "$BK/data" "$BK/overleaf"
FILES=(_calibration_override.mat RW_shock.csv data/MS_params.csv data/MS_sigma_us_prob.csv data/initguess.mat)
TEX=(Mod_Moments_cd.tex Mod_CIP_Moments_cd.tex Mod_SwitchDev_Moments_cd.tex)
for f in "${FILES[@]}"; do [ -f "$f" ] && cp "$f" "$BK/$f"; done
for f in "${TEX[@]}"; do [ -f "$OVERLEAF/$f" ] && cp "$OVERLEAF/$f" "$BK/overleaf/$f"; done
restore() { for f in "${FILES[@]}"; do if [ -f "$BK/$f" ]; then cp "$BK/$f" "$f"; fi; done
            for f in "${TEX[@]}"; do if [ -f "$BK/overleaf/$f" ]; then cp "$BK/overleaf/$f" "$OVERLEAF/$f"; fi; done
            echo "[restored]"; }
trap restore EXIT
cp _calibration_override_cbase.mat _calibration_override.mat
cp RW_shock_cbase.csv RW_shock.csv
cp data/MS_params_ms_end2021.csv data/MS_params.csv
cp data/MS_sigma_us_prob_ms_end2021.csv data/MS_sigma_us_prob.csv
cp data/initguess_cbase.mat data/initguess.mat
"$MATLAB" -batch "matching_type=1; save('temp_matching_type.mat','matching_type'); main_LFX; quit force" > end2021_lfx_log.txt 2>&1
for f in "${TEX[@]}"; do cp "$OVERLEAF/$f" "$OVERLEAF/${f%.tex}_cbase_end2021.tex"; done
echo "=== tables ==="
cat "$OVERLEAF/Mod_Moments_cd_cbase_end2021.tex"
