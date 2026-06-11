"""
Two diagnostic questions about R^f(mu, sigma; eta):

Q1: Can eta keep sigma flat when mu rises?
    - Pick a reference (mu0, sigma0, eta0) -> R^f_target.
    - For each mu, find the (sigma) needed to maintain R^f_target under each eta.
    - Plot sigma(mu) along the R^f = R^f_target level curve, for various eta.
    - If a curve is FLAT in mu, that eta perfectly absorbs the mu shift.

Q2: Holding mu fixed, how does sigma respond to changes in R^f, for various eta?
    - Plot sigma(R^f) at fixed mu, faceted by eta.
    - The slope = dsigma/dR^f -- the "amplification" the filter applies.

Outputs:
  - ted_eta_sigma_along_levelset.pdf
  - ted_sigma_of_Rf_fixed_mu.pdf
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

# Re-import the translated machinery
import importlib.util, pathlib
spec = importlib.util.spec_from_file_location(
    "core", pathlib.Path(__file__).with_name("diag_ted_surface.py")
)
core = importlib.util.module_from_spec(spec); spec.loader.exec_module(core)

# Convenience: scalar R^f
def Rf(mu, sigma, eta,
       ploss=core.PLOSS, lam=core.LAMBDA, iota=core.IOTA, varrho=core.VARRHO):
    mu_a = np.array([mu], dtype=float); sig_a = np.array([sigma], dtype=float)
    val = core.ted_surface(mu_a, sig_a, eta, ploss, lam, iota, varrho)[0]
    return float(val)

# Solve for sigma s.t. Rf(mu, sigma, eta) = R_target  (bracket on sigma > 0)
def sigma_for_Rf(mu, eta, R_target, sig_lo=1e-3, sig_hi=10.0, npts=300):
    """Solve Rf(mu, sigma, eta) = R_target for sigma. Skip NaN regions."""
    f = lambda s: Rf(mu, s, eta) - R_target
    sigs = np.linspace(sig_lo, sig_hi, npts)
    vals = np.array([f(s) for s in sigs])
    finite = np.isfinite(vals)
    if not finite.any():
        return np.nan
    # Find first interior bracket where both endpoints are finite and signs differ
    for i in range(len(vals) - 1):
        if finite[i] and finite[i+1] and (vals[i] * vals[i+1] < 0):
            try:
                return brentq(f, sigs[i], sigs[i+1])
            except Exception:
                continue
    return np.nan

# ----------------------------------------------------------------------------
# Q1: Level-set sigma(mu) under different eta for fixed R^f_target
# ----------------------------------------------------------------------------
# Reference operating point: pick a realistic post-2008/2020 sample target.
# Pre-2020 data avg R^f ≈ 240 bps (typical TED/CIP level pre-COVID).
# Choose this as the target and ask: as mu varies, how must sigma move per eta?
R_target_list = [150.0, 250.0, 350.0]      # bps/yr
mu_grid       = np.linspace(0.02, 0.5, 50)
eta_list      = [0.3, 0.5, 0.7]

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)

for k, R_target in enumerate(R_target_list):
    ax = axes[k]
    for eta in eta_list:
        sig_curve = np.array([sigma_for_Rf(mu, eta, R_target) for mu in mu_grid])
        ax.plot(mu_grid, sig_curve, lw=2, label=rf"$\eta = {eta}$")
    ax.set_xlabel(r"$\mu$ (reserve ratio)")
    if k == 0:
        ax.set_ylabel(r"required $\sigma$ to hold $R^f$ constant")
    ax.set_title(rf"$R^f$ target = {R_target:.0f} bps/yr")
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 2.5)
    ax.legend(loc='upper left', fontsize=10)

fig.suptitle(
    r"Q1: As $\mu$ rises, how much must $\sigma$ rise to hold $R^f$ fixed? "
    r"Flatter curve ⇒ $\eta$ better absorbs $\mu$-shift.",
    fontsize=12
)
plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig("ted_eta_sigma_along_levelset.pdf", bbox_inches='tight')
print("Wrote ted_eta_sigma_along_levelset.pdf")

# Print summary: at each fixed eta, what is dsigma/dmu (the slope)?
print("\n=== Q1 numerical slopes: dsigma/dmu along R^f level-set ===")
print(f"{'R^f tgt':>9} {'eta':>5} {'sigma@mu=0.05':>15} {'sigma@mu=0.20':>15} "
      f"{'(sigma_hi-sigma_lo)/(mu_hi-mu_lo)':>34}")
for R_target in R_target_list:
    for eta in eta_list:
        s_lo = sigma_for_Rf(0.05, eta, R_target)
        s_hi = sigma_for_Rf(0.20, eta, R_target)
        slope = (s_hi - s_lo) / (0.20 - 0.05) if (np.isfinite(s_lo) and np.isfinite(s_hi)) else np.nan
        print(f"{R_target:9.0f} {eta:5.1f} {s_lo:15.3f} {s_hi:15.3f} {slope:34.2f}")

# ----------------------------------------------------------------------------
# Q2: sigma(R^f) at fixed mu, for various eta
# ----------------------------------------------------------------------------
mu_fixed_list = [0.05, 0.15, 0.30]    # three reserve ratios spanning the data
fig2, axes2 = plt.subplots(1, 3, figsize=(15, 4.5))

# For each mu, compute the curve sigma(R^f) by solving inversion at many target R^f
R_grid = np.linspace(20, 500, 60)     # bps/yr

for k, mu0 in enumerate(mu_fixed_list):
    ax = axes2[k]
    for eta in eta_list:
        sig_inv = np.array([sigma_for_Rf(mu0, eta, R) for R in R_grid])
        ax.plot(R_grid, sig_inv, lw=2, label=rf"$\eta = {eta}$")
    ax.set_xlabel(r"$R^f$ (bps/yr)")
    if k == 0:
        ax.set_ylabel(r"implied $\sigma$")
    ax.set_title(rf"$\mu = {mu0}$")
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 2.5)
    ax.legend(loc='upper left', fontsize=10)

fig2.suptitle(
    r"Q2: At fixed $\mu$, what $\sigma$ does the model imply for each observed $R^f$? "
    r"Steeper curve ⇒ filter amplifies $R^f$-noise into $\sigma$-noise.",
    fontsize=12
)
plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig("ted_sigma_of_Rf_fixed_mu.pdf", bbox_inches='tight')
print("Wrote ted_sigma_of_Rf_fixed_mu.pdf")

# Print Q2 numerical slopes
print("\n=== Q2 numerical slopes: dsigma/dR^f at fixed mu ===")
print(f"{'mu':>5} {'eta':>5} {'R^f=100->150':>15} {'R^f=200->250':>15} {'R^f=300->350':>15}")
for mu0 in mu_fixed_list:
    for eta in eta_list:
        out = []
        for R0 in (100, 200, 300):
            s0 = sigma_for_Rf(mu0, eta, R0)
            s1 = sigma_for_Rf(mu0, eta, R0+50)
            if np.isfinite(s0) and np.isfinite(s1):
                out.append(f"{(s1-s0)/50:15.4f}")
            else:
                out.append(f"{'NaN':>15}")
        print(f"{mu0:5.2f} {eta:5.1f} {out[0]} {out[1]} {out[2]}")
