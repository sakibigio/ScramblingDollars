#!/usr/bin/env python3
"""
Compare the digitized FEDS-Note HQLA/total-assets series (by bank group)
against the empirical model liquidity ratio exp(mu) used by the filter.

  FEDS-Note groups : hqla_assets_fedsnote_monthly_linear_extended.csv
                     (Standard_LCR / Modified_LCR / Non_LCR banks, percent)
  Model ratio      : exp(mu_us), exp(mu_eu) from data/LFX_data.mat (fraction -> percent)

Output: data/hqla_vs_mu_linear.png
"""
import os
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import date

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")

# --- FEDS-Note HQLA/assets (linear extended), skip comment lines ----------
csv_path = os.path.join(DATA, "hqla_assets_fedsnote_monthly_pchip_extended.csv")
dates_h, std, mod, non = [], [], [], []
with open(csv_path) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("Date"):
            continue
        d, a, b, c = line.split(",")
        y, m = d.split("-")
        dates_h.append(date(int(y), int(m), 1))
        std.append(float(a)); mod.append(float(b)); non.append(float(c))
dates_h = np.array(dates_h)
std, mod, non = map(np.array, (std, mod, non))

# --- Empirical liquidity ratio exp(mu) from the filter input --------------
mat = scipy.io.loadmat(os.path.join(DATA, "LFX_data.mat"))
mu_us = np.asarray(mat["mu_us"]).ravel()
mu_eu = np.asarray(mat["mu_eu"]).ravel()
T = len(mu_us)
dates_m = np.array([date(2001 + (k)//12, (k) % 12 + 1, 1) for k in range(T)])
exp_mu_us = np.exp(mu_us) * 100.0   # fraction -> percent
exp_mu_eu = np.exp(mu_eu) * 100.0

# --- Plot -----------------------------------------------------------------
plt.rcParams.update({"font.family": "serif", "font.size": 13})
fig, ax = plt.subplots(figsize=(10, 6))

# FEDS-Note bank groups (HQLA / total assets)
ax.plot(dates_h, std, color="#1f77b4", lw=2.0, label="FEDS-Note: Standard-LCR banks")
ax.plot(dates_h, mod, color="#2ca02c", lw=2.0, label="FEDS-Note: Modified-LCR banks")
ax.plot(dates_h, non, color="#7f7f7f", lw=2.0, label="FEDS-Note: Non-LCR banks")

# Model liquidity ratio exp(mu)
ax.plot(dates_m, exp_mu_us, color="#d62728", lw=2.4, label=r"Model $e^{\mu_{US}}$ (US liquidity ratio)")
ax.plot(dates_m, exp_mu_eu, color="#d62728", lw=1.4, ls="--", alpha=0.7,
        label=r"Model $e^{\mu_{EU}}$ (EU liquidity ratio)")

# Shade the digitized (non-padded) FEDS-Note window 2003Q1-2018Q3
ax.axvspan(date(2003, 1, 1), date(2018, 9, 1), color="0.92", zorder=0,
           label="FEDS-Note digitized window")

ax.set_ylabel("Liquidity ratio  (HQLA / assets,  percent)")
ax.set_xlabel("")
ax.set_title("FEDS-Note HQLA/assets (pchip) vs. model liquidity ratio $e^{\\mu}$")
ax.xaxis.set_major_locator(mdates.YearLocator(2))
ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
ax.set_xlim(date(2001, 1, 1), date(2026, 1, 1))
ax.grid(True, alpha=0.25)
ax.legend(loc="upper left", fontsize=10, framealpha=0.9)
fig.tight_layout()

out = os.path.join(DATA, "hqla_vs_mu_pchip.png")
fig.savefig(out, dpi=150)
print("wrote", out)
print("exp(mu_us) %% range [%.1f, %.1f]; FEDS Standard range [%.1f, %.1f]"
      % (exp_mu_us.min(), exp_mu_us.max(), std.min(), std.max()))
