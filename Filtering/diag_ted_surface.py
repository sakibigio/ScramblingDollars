"""
Surface plots of R^f / TED(mu, sigma) under Cobb-Douglas matching, for varying eta.

Translates the MATLAB chain:
  distribution_aggregates(mu, ploss, sigma, varrho) -> theta
  bartheta(theta, lambda, matching_type=1) -> bar_theta
  Chi_p(theta, iota, lambda, eta, matching_type=1) -> chi_p
  Psi_p(theta, lambda, matching_type=1)       -> psi_p
  TED = chi_p / psi_p

Output: ted_surface_eta.pdf and ted_slice_eta.pdf
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# ---------------------------------------------------------------------------
# Calibration (Cobb-Douglas — matches params.m current values)
# ---------------------------------------------------------------------------
PLOSS    = 0.5
LAMBDA   = 1.2778            # CD matching efficiency
IOTA_ANN = 0.0625
IOTA     = IOTA_ANN / 12     # monthly (DW-IOR spread)
VARRHO   = 0
MATCHING = 1                 # 1 = Cobb-Douglas

# Display scaling (annualized basis points, as used in main_filter.m line 257)
ABS_SCALE = 1e4 * 12

# ---------------------------------------------------------------------------
# Translated helpers (vectorized over arrays)
# ---------------------------------------------------------------------------
def distribution_aggregates(mu, ploss, sigma, varrho=0):
    """Returns (Smin, Spl, theta) for asymmetric Laplace shock distribution."""
    m_eff = mu - varrho
    expo_lo = sigma * np.exp(-ploss     / sigma * m_eff)   # for m_eff >= 0
    expo_hi = sigma * np.exp((1 - ploss) / sigma * m_eff)  # for m_eff <  0
    idx_pos = (m_eff >= 0)
    Smin = np.where(idx_pos, expo_lo,        -m_eff + expo_hi)
    Spl  = np.where(idx_pos, m_eff + expo_lo, expo_hi)
    theta = Smin / Spl
    return Smin, Spl, theta

def bartheta_CD(theta, lam):
    """Cobb-Douglas terminal tightness."""
    bar = np.zeros_like(theta)
    eps = 1e-12

    # theta < 1 branch
    mask_lo = theta < 1 - eps
    if mask_lo.any():
        th = theta[mask_lo]
        sqrt_th = np.sqrt(th)
        alpha = (1 + sqrt_th) / (1 - sqrt_th)
        T_star = np.log(np.abs(alpha)) / lam
        T_max = np.minimum(T_star, 1.0)
        exp_T = alpha * np.exp(-lam * T_max)
        val = ((exp_T - 1) / (exp_T + 1))**2
        val = np.where(T_star < 1, 0.0, val)   # early stop
        bar[mask_lo] = val

    # theta > 1 branch
    mask_hi = theta > 1 + eps
    if mask_hi.any():
        th = theta[mask_hi]
        sqrt_th = np.sqrt(th)
        # alpha for theta > 1: numerator/denominator both real if formula applied carefully
        alpha = (1 + sqrt_th) / (1 - sqrt_th)   # negative
        T_star = np.log(np.abs(alpha)) / lam
        T_max = np.minimum(T_star, 1.0)
        exp_T = alpha * np.exp(-lam * T_max)
        val = ((exp_T - 1) / (exp_T + 1))**2
        val = np.where(T_star < 1, np.inf, val)
        bar[mask_hi] = val

    mask_one = np.abs(theta - 1) <= eps
    bar[mask_one] = 1.0
    return bar

def Chi_p_CD(theta, iota, lam, eta):
    """Liquidity yield chi_p for Cobb-Douglas matching."""
    barth = bartheta_CD(theta, lam)
    chip = np.zeros_like(theta)
    idx_gen = np.abs(barth - 1) > 1e-10
    idx_one = ~idx_gen
    chip[idx_gen] = iota * (barth[idx_gen]
                             - barth[idx_gen]**eta * theta[idx_gen]**(1 - eta)) \
                          / (barth[idx_gen] - 1)
    chip[idx_one] = iota * (1 - eta)
    return chip

def Psi_p_CD(theta, lam):
    """Matching probability Psi_p for Cobb-Douglas matching."""
    psi = np.zeros_like(theta)
    # piecewise per theta vs 1
    for mask in (theta < 1, theta > 1):
        if not mask.any():
            continue
        th = theta[mask]
        sqrt_th = np.sqrt(th)
        alpha = (1 + sqrt_th) / (1 - sqrt_th)
        T_star = np.log(np.abs(alpha)) / lam
        T_max = np.minimum(T_star, 1.0)
        val = 1 - np.exp(-lam * T_max) * \
                  ((alpha + np.exp(lam * T_max)) / (alpha + 1))**2
        psi[mask] = val
    psi[theta == 1] = 1 - np.exp(-lam)
    return psi

def ted_surface(MU, SIG, eta,
                ploss=PLOSS, lam=LAMBDA, iota=IOTA, varrho=VARRHO):
    """Vectorized TED on grids MU, SIG. Returns (TED * ABS_SCALE)."""
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        _, _, theta = distribution_aggregates(MU, ploss, SIG, varrho)
        chip = Chi_p_CD(theta, iota, lam, eta)
        psip = Psi_p_CD(theta, lam)
        ted = np.where(psip > 0, chip / psip, np.nan)
    return ted * ABS_SCALE

def safe_q(arr, q=99):
    """Quantile robust to nans/infs."""
    finite = arr[np.isfinite(arr)]
    return float(np.percentile(finite, q)) if finite.size else 1.0

# ---------------------------------------------------------------------------
# Build grid
# ---------------------------------------------------------------------------
# mu range: reserve ratio is in [0, 1]. Most action happens in [0, 0.5].
mu_vec = np.linspace(0.0, 1.0, 100)
# sigma range: covers data range plus headroom
sigma_vec = np.linspace(0.05, 1.5, 80)
MU, SIG = np.meshgrid(mu_vec, sigma_vec, indexing='xy')

eta_values = [0.3, 0.5, 0.7, 0.9]

# ---------------------------------------------------------------------------
# Plot 1: 2×2 surface plots — TED(mu, sigma) for each eta
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(13, 10))
fig.suptitle(
    r"$R^f$ (TED, annualized bps) as a function of $\mu$ and $\sigma$ — Cobb–Douglas, "
    f"$\\lambda$={LAMBDA}, $\\iota$={IOTA_ANN}/yr, $p$={PLOSS}",
    fontsize=12, y=0.995
)

# Determine common vertical scale across eta panels (use 99% quantile to avoid one-cell blowups)
all_max = max(safe_q(ted_surface(MU, SIG, e), 99) for e in eta_values)

for k, eta in enumerate(eta_values):
    ax = fig.add_subplot(2, 2, k+1, projection='3d')
    TED = ted_surface(MU, SIG, eta)
    TED_clean = np.where(np.isfinite(TED), np.minimum(TED, all_max), np.nan)
    surf = ax.plot_surface(MU, SIG, TED_clean,
                           cmap=cm.viridis, edgecolor='none',
                           antialiased=True, vmin=0, vmax=all_max)
    ax.set_xlabel(r"$\mu$ (reserve ratio)")
    ax.set_ylabel(r"$\sigma$ (vol)")
    ax.set_zlabel(r"$R^f$ (bps/yr)")
    ax.set_title(rf"$\eta$ = {eta}", fontsize=12)
    ax.view_init(elev=22, azim=-130)
    ax.set_zlim(0, all_max)

cbar = fig.colorbar(surf, ax=fig.axes, shrink=0.45, aspect=15, pad=0.06,
                    label=r"$R^f$ (bps/yr, clipped at 99-pct)")
plt.savefig("ted_surface_eta.pdf", bbox_inches='tight')
print(f"Wrote ted_surface_eta.pdf (z clipped at {all_max:.0f} bps/yr)")

# ---------------------------------------------------------------------------
# Plot 2: 2D slices — TED vs sigma for fixed mu, faceted by eta
# Shows the σ-sensitivity directly. Curves grouped by mu value.
# ---------------------------------------------------------------------------
fig2, axes = plt.subplots(1, 4, figsize=(18, 4.5), sharey=True)
mu_slices = [0.02, 0.05, 0.1, 0.3, 0.6]
sigma_fine = np.linspace(0.05, 1.5, 300)

for k, eta in enumerate(eta_values):
    ax = axes[k]
    for mu0 in mu_slices:
        MU0 = np.full_like(sigma_fine, mu0)
        TED = ted_surface(MU0, sigma_fine, eta)
        ax.plot(sigma_fine, TED, label=rf"$\mu={mu0}$", lw=2)
    ax.set_xlabel(r"$\sigma$")
    if k == 0:
        ax.set_ylabel(r"$R^f$ (bps/yr)")
    ax.set_title(rf"$\eta$ = {eta}")
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 800)

axes[-1].legend(loc="upper right", fontsize=10)
fig2.suptitle(
    rf"$R^f(\mu, \sigma)$ slices — sensitivity to $\sigma$ at fixed $\mu$, "
    rf"varying $\eta$. CD, $\lambda$={LAMBDA}, $\iota$={IOTA_ANN}/yr, $p$={PLOSS}",
    fontsize=12
)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig("ted_slice_eta.pdf", bbox_inches='tight')
print("Wrote ted_slice_eta.pdf")

# ---------------------------------------------------------------------------
# Plot 3: sensitivity dR^f / dsigma vs sigma (numerical derivative)
# Most explicit picture of "where in parameter space is the sensitivity high"
# ---------------------------------------------------------------------------
fig3, axes3 = plt.subplots(1, 4, figsize=(18, 4.5), sharey=True)
sigma_grid = np.linspace(0.05, 1.5, 400)
dsig = sigma_grid[1] - sigma_grid[0]

for k, eta in enumerate(eta_values):
    ax = axes3[k]
    for mu0 in mu_slices:
        MU0 = np.full_like(sigma_grid, mu0)
        TED = ted_surface(MU0, sigma_grid, eta)
        # central difference
        dTED = np.gradient(TED, dsig)
        ax.plot(sigma_grid, dTED, label=rf"$\mu={mu0}$", lw=2)
    ax.set_xlabel(r"$\sigma$")
    if k == 0:
        ax.set_ylabel(r"$dR^f/d\sigma$ (bps/yr per unit $\sigma$)")
    ax.set_title(rf"$\eta$ = {eta}")
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='k', lw=0.5)

axes3[-1].legend(loc="upper right", fontsize=10)
fig3.suptitle(
    r"Sensitivity $dR^f/d\sigma$ — where in $(\mu,\sigma)$ space does $\sigma$ "
    r"have the biggest effect on $R^f$?",
    fontsize=12
)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig("ted_dsigma_eta.pdf", bbox_inches='tight')
print("Wrote ted_dsigma_eta.pdf")

# ---------------------------------------------------------------------------
# Print useful summary numbers
# ---------------------------------------------------------------------------
print("\n=== Summary at the data-relevant operating point ===")
print("Recent σ_us drift in filter: ~0.20 (post-2020 plateau)")
print("Pre-2020 σ_us range:        0.14 to 3.04 (most mass 0.20–0.55)")
print()
print(f"{'eta':>5} {'mu':>6} {'sigma':>7} {'TED bps':>9} {'dTED/dsig':>10}")
for eta in eta_values:
    for mu0 in [0.02, 0.05, 0.10, 0.30]:
        for s0 in [0.20, 0.30, 0.50]:
            t1 = ted_surface(np.array([mu0]), np.array([s0]),       eta)[0]
            t2 = ted_surface(np.array([mu0]), np.array([s0+1e-3]),  eta)[0]
            print(f"{eta:5.1f} {mu0:6.2f} {s0:7.3f} {t1:9.1f} {(t2-t1)/1e-3:10.1f}")
