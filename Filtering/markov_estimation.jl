# Markov Switching Estimation for Sigma Shocks
# Part of: Scrambling for Dollars pipeline
#
# Reads: RW_shock.csv (from main_filter.m)
# Writes:
#   data/MS_sigma_us_prob.csv          - filtered regime probabilities (cols: prob_nor, prob_scr)
#   data/MS_sigma_us_counterfactuals.csv - counterfactual sigma path
#   data/MS_sigma_us_simul.csv         - simulated paths
#   data/MS_sigma_eu_*.csv             - same for EU
#   data/MS_params.csv                 - estimated parameters for setup_markov.m
#                                        r1 = normal, r2 = scrambling (consistently)
#
# Split-sample design (see EST_END_IDX below):
#   - Parameters are ESTIMATED on the pre-2020 subsample only, because the
#     post-2020 filter output is unusually flat and distorts regime identification.
#   - Filtered probabilities are REPORTED over the FULL sample (including
#     post-2020) using a Hamilton forward filter with frozen pre-2020 parameters.
#   - Counterfactuals and simulations use the pre-2020 parameters.

using MarSwitching
using CSV
using DataFrames
using Plots
using Random
using Distributions
using StatsBase
using Statistics
using LinearAlgebra   # eigen() for Hamilton-filter stationary init

## =========================================================================
##  CONFIG: split estimation vs. probability-inference samples
## =========================================================================
## RW_shock.csv is monthly starting 2001-01. After dropmissing on the lag
## variables, df row 1 = 2001-02 and df row 227 = 2019-12.
##
## The post-2020 portion of the filtered σ_us series is unusually flat
## (σ_us drifts from ~0.27 to ~0.20 over 2021–2025 with std ≈ 0.06).
## Including it in estimation forces a near-unit-root ρ_r1 and collapses
## the μ-separation between normal and scrambling regimes. We therefore:
##
##   * ESTIMATE parameters on df[1:EST_END_IDX, :]  (pre-2020 only)
##   * REPORT filtered probabilities for the FULL sample, using a Hamilton
##     forward filter with the frozen pre-2020 parameters.
##
## Set EST_END_IDX = nrow(df) to revert to full-sample estimation.
## Env override: MS_EST_END=<row>  (MS_EST_END=0 -> full sample). Default 227
## keeps the committed pre-2020 truncation.
const EST_END_IDX = let v = get(ENV, "MS_EST_END", "")
    isempty(v) ? 227 : parse(Int, v)   # df-row of last estimation observation (227 = Dec 2019)
end

## Stationarity screen for the multistart (env MS_RHO_MAX, default Inf = off).
## When set (e.g. 0.99), seeds whose fitted regime AR coefficients include
## max_k rho_k >= MS_RHO_MAX are rejected as degenerate, so the 21-seed search
## returns the best STATIONARY local optimum instead of the explosive global
## one. Rationale: regime long-run means (intercept/(1-rho)) are undefined for
## nonstationary regimes, so such fits are uninterpretable, not just inferior.
const RHO_MAX = let v = get(ENV, "MS_RHO_MAX", "")
    isempty(v) ? Inf : parse(Float64, v)
end

## μ as additional switching regressor (orthogonalization)
## When true, the AR(1) becomes:
##   log σ_t = c_k + ρ_k · log σ_{t-1} + γ_k · μ_t + ε_{k,t}
## which is equivalent to the user's parameterization
##   log σ_t = c_k + α_k·(1-ρ_k)·μ_t + ρ_k · log σ_{t-1} + ε_{k,t}
## with α_k = γ_k / (1 - ρ_k) being the long-run μ-loading per state.
## Requires RW_shock.csv to contain a `mu_us` column (and `mu_eu` for the EU model).
const USE_MU_REGRESSOR = false

## When USE_MU_TREND is true, the μ-regressor is replaced by its Hodrick-Prescott
## trend component (HP_lambda below). Idea: capture the secular drift in reserves
## as the long-run regressor, leaving the cyclical μ variation in the residuals.
## NOTE: kept false in the locked-in version (raw μ produced the published moments).
const USE_MU_TREND  = false
const HP_LAMBDA_MO  = 129600.0   # Ravn-Uhlig for monthly data

## MU_SHARED_GAMMA: when true, μ enters as a NON-SWITCHING regressor — one
## common γ across both regimes. AR(1) becomes:
##   log σ_t = c_k + ρ_k · log σ_{t-1} + γ · μ_t + ε_{k,t}
## (γ shared; c_k and ρ_k still regime-dependent). Parsimonious version.
## When false, γ is regime-specific (the original spec).
const MU_SHARED_GAMMA = true

## USE_LCR_REGRESSOR: US-only alternative to the μ-regressor. When true, the US
## model uses the empirical LCR/HQLA series (`lcr_us` column in RW_shock.csv) as
## the orthogonalization regressor IN PLACE OF μ — soaking up the regulatory
## liquidity trend that drives the σ–μ coupling, without using μ itself. Routed
## through the same machinery as the μ-regressor (non-switching when
## MU_SHARED_GAMMA=true). EU is unaffected (no EU LCR series). Mutually exclusive
## with USE_MU_REGRESSOR for the US model.
const USE_LCR_REGRESSOR = false

## When true, the LCR regressor enters as log(LCR/100) — a LOG liquidity ratio,
## matching the form μ takes under USE_MU_REGRESSOR (df.mu_us is a log ratio).
## When false, the LCR enters in raw level (percent); for a linear regressor the
## level scale is absorbed by γ, so percent vs ratio is immaterial — only
## level-vs-log changes the fit.
const LCR_AS_LOG_RATIO = true

# Two-sided HP filter via direct linear solve.
# Returns trend τ minimizing Σ (y_t - τ_t)² + λ · Σ (τ_{t+1} - 2τ_t + τ_{t-1})²
function hp_trend(y::AbstractVector, lam::Real)
    T = length(y)
    # Second-difference matrix D: (T-2)×T
    D = zeros(T-2, T)
    for i in 1:T-2
        D[i, i]   =  1.0
        D[i, i+1] = -2.0
        D[i, i+2] =  1.0
    end
    A = I + lam * (D' * D)
    return A \ y
end

## Load data (relative path - run from Filtering folder)
df = CSV.read("RW_shock.csv", DataFrame, missingstring = "NA")

# Create lagged variables
df.lag_sigma_us = [missing; df.sigma_us[1:end-1]]
df.lag_sigma_eu = [missing; df.sigma_eu[1:end-1]]
dropmissing!(df, :lag_sigma_eu)
dropmissing!(df, :lag_sigma_us)

# Log transform
df.lsigma_us = log.(df.sigma_us)
df.lag_lsigma_us = log.(df.lag_sigma_us)
df.lsigma_eu = log.(df.sigma_eu)
df.lag_lsigma_eu = log.(df.lag_sigma_eu)

# Plot raw data
p1 = plot(df.sigma_us, label="σ_us", title="Filtered Sigma Shocks")
plot!(p1, df.sigma_eu, label="σ_eu")
savefig(p1, "data/sigma_raw.png")

## =========================================================================
##  Robust estimation: retry with different seeds, reject degenerate runs
## =========================================================================

function is_degenerate(model; min_sigma=0.01, min_prob_balance=0.02)
    if minimum(model.σ) < min_sigma
        return true, "min σ = $(minimum(model.σ)) < $min_sigma"
    end
    # Stationarity screen (off unless MS_RHO_MAX is set). The AR(1) coefficient
    # is the first switching regressor after the intercept: model.β[k][2].
    if isfinite(RHO_MAX) && all(length(b) >= 2 for b in model.β)
        rho_k = [model.β[k][2] for k in eachindex(model.β)]
        if maximum(rho_k) >= RHO_MAX
            return true, "max ρ = $(round(maximum(rho_k), digits=4)) ≥ $RHO_MAX (nonstationary regime)"
        end
    end
    try
        fp = filtered_probs(model)
        fp_clean = fp[.!any(isnan.(fp), dims=2)[:], :]
        if size(fp_clean, 1) < size(fp, 1) * 0.5
            return true, "$(size(fp,1) - size(fp_clean,1)) NaN rows in filtered probs"
        end
        mean_p1 = mean(fp_clean[:,1])
        if mean_p1 < min_prob_balance || mean_p1 > (1 - min_prob_balance)
            return true, "mean(prob_state1) = $(round(mean_p1, digits=4)) — one state never visited"
        end
    catch e
        return true, "filtered_probs error: $e"
    end
    return false, ""
end

function robust_msmodel(y, n_states, exog_sv; label="", max_retries=20, min_sigma=0.01, exog_ns=nothing)
    best_model = nothing
    best_ll = -Inf

    for i in 0:max_retries
        if i > 0
            Random.seed!(i)
        end
        try
            m = if exog_ns === nothing
                MSModel(y, n_states, exog_switching_vars = exog_sv)
            else
                MSModel(y, n_states, exog_switching_vars = exog_sv, exog_vars = exog_ns)
            end
            ll = m.Likelihood

            if isnan(ll)
                println("  $label seed $i: rejected (loglik=NaN)")
                continue
            end

            degen, reason = is_degenerate(m; min_sigma=min_sigma)
            if degen
                println("  $label seed $i: degenerate ($reason), loglik=$ll")
                continue
            end

            println("  $label seed $i: loglik=$(round(ll, digits=2)), σ=$(round.(m.σ, digits=4))")
            if ll > best_ll
                best_ll = ll
                best_model = m
            end
        catch e
            println("  $label seed $i: error — $e")
        end
    end

    if best_model === nothing
        println("  $label: all $(max_retries+1) attempts degenerate or failed")
        return nothing
    end

    println("\n  $label BEST: loglik=$(round(best_ll, digits=2))")
    println("    σ = $(round.(best_model.σ, digits=4))")
    println("    β₁ = $(round.(best_model.β[1], digits=4))")
    println("    β₂ = $(round.(best_model.β[2], digits=4))")
    fp = filtered_probs(best_model)
    fp_clean = fp[.!any(isnan.(fp), dims=2)[:], :]
    println("    mean(prob_state1) = $(round(mean(fp_clean[:,1]), digits=3))")
    println("    non-NaN filtered prob rows: $(size(fp_clean,1)) / $(size(fp,1))")

    return best_model
end

## =========================================================================
##  Regime labelling helpers
## =========================================================================
##
## Convention used throughout this script and downstream:
##   nor (r1) = "normal"     = the MORE frequent regime (higher mean filtered prob)
##   scr (r2) = "scrambling" = the LESS frequent regime
##
## MarSwitching returns model.P as a COLUMN-stochastic matrix:
##   P_raw[i, j] = Pr(state_t = i | state_{t-1} = j)
## We re-emit P as ROW-stochastic, with rows/cols ordered (nor, scr):
##   P_out[r, c] = Pr(state_t = state_c | state_{t-1} = state_r)
## Using the same (idx_nor, idx_scr) ordering for β, σ AND P prevents any
## mismatch between AR(1) coefficients and transition probabilities.

function regime_indices(model)
    fp = filtered_probs(model)
    if mean(fp[:,1]) >= mean(fp[:,2])
        return 1, 2   # idx_nor, idx_scr
    else
        return 2, 1
    end
end

function reorder_P(P_raw, idx_nor, idx_scr)
    # Convert column-stochastic P_raw to row-stochastic with rows (nor, scr)
    return [P_raw[idx_nor, idx_nor]  P_raw[idx_scr, idx_nor];
            P_raw[idx_nor, idx_scr]  P_raw[idx_scr, idx_scr]]
end

## =========================================================================
##  Hamilton forward filter with FROZEN parameters
## =========================================================================
## Given a fitted MSModel (parameters estimated on any subsample) and an
## arbitrary AR(1) observation series (y_t, y_{t-1}), runs Hamilton's forward
## filter to produce filtered regime probabilities in (nor, scr) ordering.
##
## Used here to back out FULL-SAMPLE probabilities from PRE-2020 estimates.
function hamilton_filter_ar1(y::AbstractVector, X_sw::AbstractVecOrMat,
                              model, idx_nor::Int, idx_scr::Int;
                              X_ns::Union{AbstractMatrix,Nothing}=nothing)
    # MarSwitching stores ALL coefficients (switching + non-switching) in
    # model.β[k] in order [intercept, switching_vars..., non_switching_vars...].
    # For non-switching coefs, β[k][last] is identical across k.
    # We concatenate X_sw and X_ns and just use β[k][1] + β[k][2:end]·X[t,:].
    T = length(y)
    if X_sw isa AbstractVector
        X_sw_mat = reshape(X_sw, T, 1)
    else
        X_sw_mat = X_sw
    end
    @assert size(X_sw_mat, 1) == T "X_sw must have T rows"
    Xmat = X_ns === nothing ? X_sw_mat : hcat(X_sw_mat, X_ns)

    β_nor = model.β[idx_nor];  β_scr = model.β[idx_scr]
    σ_nor = model.σ[idx_nor];  σ_scr = model.σ[idx_scr]
    P_rs  = reorder_P(model.P, idx_nor, idx_scr)

    ev = eigen(P_rs')
    j  = argmin(abs.(ev.values .- 1.0))
    π0 = real.(ev.vectors[:, j]);  π0 = π0 ./ sum(π0)

    prob_nor = zeros(T)
    prob_scr = zeros(T)
    ξ = π0
    for t in 1:T
        ξ_pred = P_rs' * ξ
        m_nor_t = β_nor[1] + sum(β_nor[i+1] * Xmat[t, i] for i in 1:size(Xmat,2))
        m_scr_t = β_scr[1] + sum(β_scr[i+1] * Xmat[t, i] for i in 1:size(Xmat,2))
        f_nor = pdf(Normal(m_nor_t, σ_nor), y[t])
        f_scr = pdf(Normal(m_scr_t, σ_scr), y[t])

        num   = (ξ_pred[1] * f_nor, ξ_pred[2] * f_scr)
        denom = num[1] + num[2]
        if denom == 0 || !isfinite(denom)
            ξ = ξ_pred
        else
            ξ = [num[1] / denom, num[2] / denom]
        end

        prob_nor[t] = ξ[1]
        prob_scr[t] = ξ[2]
    end
    return prob_nor, prob_scr
end

## =========================================================================
##  US Markov Switching Model
## =========================================================================
println("\n=== US Markov Switching Model ===")

# Estimation subsample (pre-2020 by default; MS_EST_END=0 -> full sample).
EST_END = (EST_END_IDX <= 0 || EST_END_IDX > nrow(df)) ? nrow(df) : EST_END_IDX
df_est = df[1:EST_END, :]
println("Estimation sample: rows 1..$(EST_END) of $(nrow(df))" *
        (EST_END < nrow(df) ? "  (truncated)" : "  (FULL SAMPLE)"))
println("Probability inference (Hamilton filter): full sample, rows 1..$(nrow(df))")
# US regressor source: μ (USE_MU_REGRESSOR) or LCR (USE_LCR_REGRESSOR). LCR wins.
us_reg_on  = USE_MU_REGRESSOR || USE_LCR_REGRESSOR
us_reg_col = USE_LCR_REGRESSOR ? :lcr_us : :mu_us
reg_lbl    = USE_LCR_REGRESSOR ? "LCR" : (USE_MU_TREND ? "μ_trend" : "μ")
println("USE_MU_REGRESSOR = $(USE_MU_REGRESSOR) | USE_LCR_REGRESSOR = $(USE_LCR_REGRESSOR) (US regressor = $(reg_lbl))")

if us_reg_on && hasproperty(df, us_reg_col)
    if USE_MU_TREND && !USE_LCR_REGRESSOR
        mu_us_full  = hp_trend(df[!, us_reg_col], HP_LAMBDA_MO)
        mu_us_est   = mu_us_full[1:EST_END]
        raw_lo, raw_hi = extrema(df[!, us_reg_col]);  trd_lo, trd_hi = extrema(mu_us_full)
        println("US: HP-trend of μ used (λ=$(HP_LAMBDA_MO)). " *
                "Raw μ range=[$(round(raw_lo,digits=3)), $(round(raw_hi,digits=3))], " *
                "trend range=[$(round(trd_lo,digits=3)), $(round(trd_hi,digits=3))]")
    else
        mu_us_full  = df[!, us_reg_col]
        mu_us_est   = df_est[!, us_reg_col]
        if USE_LCR_REGRESSOR && LCR_AS_LOG_RATIO
            mu_us_full = log.(mu_us_full ./ 100)   # percent -> log liquidity ratio
            mu_us_est  = log.(mu_us_est  ./ 100)
            println("US: LCR entered as log(LCR/100) to match μ's form")
        end
    end
    if MU_SHARED_GAMMA
        # regressor enters as NON-switching (single γ across regimes)
        X_sw_us_est  = df_est.lag_lsigma_us           # only lag is switching
        X_sw_us_full = df.lag_lsigma_us
        X_ns_us_est  = reshape(collect(mu_us_est),  :, 1)   # 1-col non-switching
        X_ns_us_full = reshape(collect(mu_us_full), :, 1)
        println("US: passing lag_lsigma as switching, $(reg_lbl) as NON-switching (shared γ)")
    else
        X_sw_us_est  = hcat(df_est.lag_lsigma_us, mu_us_est)
        X_sw_us_full = hcat(df.lag_lsigma_us,     mu_us_full)
        X_ns_us_est  = nothing; X_ns_us_full = nothing
        println("US: passing 2-column switching regressors [lag_lsigma, $(reg_lbl)] (regime-specific γ)")
    end
else
    X_sw_us_est  = df_est.lag_lsigma_us
    X_sw_us_full = df.lag_lsigma_us
    X_ns_us_est  = nothing; X_ns_us_full = nothing
    if us_reg_on
        println("WARNING: US regressor requested but $(us_reg_col) column missing from RW_shock.csv -- falling back to lag-only.")
    end
end

model_us = robust_msmodel(df_est.lsigma_us, 2, X_sw_us_est,
                          label="US (pre-2020)", exog_ns=X_ns_us_est)
if model_us === nothing
    error("US: robust estimation failed — no non-degenerate solution found")
end

try
    summary_msm(model_us)
catch e
    println("Warning: Could not display US model summary (std errors singular)")
    println("  Coefficients β: ", model_us.β)
    println("  Residual σ:     ", model_us.σ)
end

# Identify regimes once, by filtered-probability frequency (on estimation sample)
idx_nor_us, idx_scr_us = regime_indices(model_us)
println("US regime mapping: normal = MarSwitching state $idx_nor_us, scrambling = state $idx_scr_us")

# Full-sample filtered probabilities via Hamilton forward filter with frozen
# pre-2020 parameters. Length = nrow(df), covering the post-2020 period too.
prob_nor_us, prob_scr_us = hamilton_filter_ar1(
    df.lsigma_us, X_sw_us_full, model_us, idx_nor_us, idx_scr_us;
    X_ns = X_ns_us_full
)
println("US filtered probabilities (full sample, Hamilton): ",
        "mean prob_scr = $(round(mean(prob_scr_us), digits=3))")

# Save filtered probabilities under semantic column names
df_probs_us = DataFrame(prob_nor = prob_nor_us, prob_scr = prob_scr_us)
CSV.write("data/MS_sigma_us_prob.csv", df_probs_us)

p2 = plot(prob_nor_us, label="Normal", title="US Regime Probabilities")
plot!(p2, prob_scr_us, label="Scrambling")
savefig(p2, "data/MS_sigma_us_probs.png")

# Extract regime-specific AR(1) parameters using the same indices
intercept_nor_us = model_us.β[idx_nor_us][1]
phi_nor_us       = model_us.β[idx_nor_us][2]
sigma_nor_us     = model_us.σ[idx_nor_us]
intercept_scr_us = model_us.β[idx_scr_us][1]
phi_scr_us       = model_us.β[idx_scr_us][2]
sigma_scr_us     = model_us.σ[idx_scr_us]

# Long-run means.
# Standard AR(1):                                     mu = c / (1 - ρ)
# With μ-regressor (β = [c, ρ, γ_μ]):                 mu = (c + γ_μ * E[μ]) / (1 - ρ)
# The added term (γ_μ * E[μ]) shifts the regime's long-run level by its μ-loading
# at the sample-average μ. setup_markov.m uses these means to build the Tauchen grid,
# so this correction matters for the discretization to be centered correctly.
if us_reg_on && hasproperty(df_est, us_reg_col) && length(model_us.β[idx_nor_us]) >= 3
    # γ is in β[k][end]: same value across regimes when MU_SHARED_GAMMA=true,
    # regime-specific when false. (regressor = μ or LCR per us_reg_col)
    mu_us_bar = mean(mu_us_est)
    gamma_nor_us = model_us.β[idx_nor_us][end]
    gamma_scr_us = model_us.β[idx_scr_us][end]
    mu_nor_us = (intercept_nor_us + gamma_nor_us * mu_us_bar) / (1 - phi_nor_us)
    mu_scr_us = (intercept_scr_us + gamma_scr_us * mu_us_bar) / (1 - phi_scr_us)
    src_lbl = reg_lbl
    src_lbl_extra = MU_SHARED_GAMMA ? " [shared γ]" : " [per-regime γ]"
    println("Long-run means (at E[$(src_lbl)_us]=$(round(mu_us_bar,digits=4)))$(src_lbl_extra): " *
            "nor=$(round(mu_nor_us,digits=4)), scr=$(round(mu_scr_us,digits=4))")
else
    mu_nor_us = intercept_nor_us / (1 - phi_nor_us)
    mu_scr_us = intercept_scr_us / (1 - phi_scr_us)
end
sigma_ss_nor_us = exp(mu_nor_us)
sigma_ss_scr_us = exp(mu_scr_us)

# Transition matrix in (nor, scr) ordering, row-stochastic
P_us_ordered = reorder_P(model_us.P, idx_nor_us, idx_scr_us)
trans_nor_us = 1 - P_us_ordered[1, 1]   # prob of leaving normal
trans_scr_us = 1 - P_us_ordered[2, 2]   # prob of leaving scrambling

# Simulate (state values follow MarSwitching's internal ordering)
simul_us           = generate_msm(model_us, 10000)
sigma_us_sim       = simul_us[1]
sigma_us_state_sim = simul_us[2]
aux_us = DataFrame(sigma_us_sim = sigma_us_sim, state_us_sim = sigma_us_state_sim)
CSV.write("data/MS_sigma_us_simul.csv", aux_us)

println("Estimated AR(1) coefficient (phi_us, normal):     ", phi_nor_us)
println("Estimated AR(1) coefficient (phi_us, scrambling): ", phi_scr_us)

# ── Standard errors via Hessian ──────────────────────────────────────────
function compute_msm_se(model, idx_nor, idx_scr, ss_nor, ss_scr, phi_nor, phi_scr)
    nans = (NaN, NaN, NaN, NaN, NaN, NaN, NaN)
    try
        raw_se   = MarSwitching.get_std_errors(model)
        V_σ, V_β = MarSwitching.vec2param_switch(raw_se, model.k,
                        model.n_β, model.n_β_ns, model.switching_var)
        se_Σ_nor   = V_σ[idx_nor];     se_Σ_scr   = V_σ[idx_scr]
        se_int_nor = V_β[idx_nor][1];  se_int_scr = V_β[idx_scr][1]
        se_ρ_nor   = V_β[idx_nor][2];  se_ρ_scr   = V_β[idx_scr][2]
        # γ (μ-loading) is the 3rd β coefficient when USE_MU_REGRESSOR=true.
        # If shared (MU_SHARED_GAMMA=true), V_β[k][3] is the same across k.
        se_gamma = length(V_β[idx_nor]) >= 3 ? V_β[idx_nor][3] : NaN
        # Delta method: SE(exp(mu)) ≈ exp(mu) * SE(intercept) / |1 - phi|
        se_ss_nor = ss_nor * se_int_nor / abs(1 - phi_nor)
        se_ss_scr = ss_scr * se_int_scr / abs(1 - phi_scr)
        println("  Standard errors computed successfully")
        return (se_ss_nor, se_ss_scr, se_ρ_nor, se_ρ_scr, se_Σ_nor, se_Σ_scr, se_gamma)
    catch e
        println("  Warning: could not compute standard errors — $e")
        return nans
    end
end

(se_sigma_ss_nor, se_sigma_ss_scr, se_rho_nor, se_rho_scr, se_Sigma_nor, se_Sigma_scr, se_gamma_us) =
    compute_msm_se(model_us, idx_nor_us, idx_scr_us,
                   sigma_ss_nor_us, sigma_ss_scr_us, phi_nor_us, phi_scr_us)

# γ value: use β[k][end] (same across k under shared-γ constraint)
gamma_us_val = (USE_MU_REGRESSOR && length(model_us.β[idx_nor_us]) >= 3) ?
                model_us.β[idx_nor_us][end] : NaN

us_params_df = DataFrame(
    param = [
        "sigma_ss_scr",    "sigma_ss_nor",
        "sigma_ss_scr_se", "sigma_ss_nor_se",
        "rho_scr",         "rho_nor",
        "rho_scr_se",      "rho_nor_se",
        "Sigma_scr",       "Sigma_nor",
        "Sigma_scr_se",    "Sigma_nor_se",
        "trans_scr",       "trans_nor",
        "gamma_mu",        "gamma_mu_se"
    ],
    value = [
        sigma_ss_scr_us,   sigma_ss_nor_us,
        se_sigma_ss_scr,   se_sigma_ss_nor,
        phi_scr_us,        phi_nor_us,
        se_rho_scr,        se_rho_nor,
        sigma_scr_us,      sigma_nor_us,
        se_Sigma_scr,      se_Sigma_nor,
        trans_scr_us,      trans_nor_us,
        gamma_us_val,      se_gamma_us
    ]
)

CSV.write("data/MS_sigma_us_params.csv", us_params_df)
println("\nUS params written to data/MS_sigma_us_params.csv:")
println(us_params_df)

# ── Generate US counterfactual (replace scrambling regimes with normal AR(1))
lsigma_us = df.lsigma_us
N = length(lsigma_us)
lsigma_us_filtered = zeros(N)
cutoff = 0.25

for t in 1:N
    if t == 1
        lsigma_us_filtered[t] = lsigma_us[t]
    elseif prob_scr_us[t-1] < cutoff && prob_scr_us[t] >= cutoff
        lsigma_us_filtered[t] = phi_nor_us * lsigma_us[t-1] + intercept_nor_us
    elseif prob_scr_us[t-1] >= cutoff && prob_scr_us[t] >= cutoff
        lsigma_us_filtered[t] = phi_nor_us * lsigma_us_filtered[t-1] + intercept_nor_us
    else
        lsigma_us_filtered[t] = lsigma_us[t]
    end
end

sigma_us_filtered = exp.(lsigma_us_filtered)

p3 = plot(exp.(lsigma_us), label="Actual", title="US Sigma: Actual vs Counterfactual")
plot!(p3, sigma_us_filtered, label="Counterfactual")
savefig(p3, "data/MS_sigma_us_counterfactual.png")

df_us_filtered = DataFrame(sigma_us_filtered = sigma_us_filtered)
CSV.write("data/MS_sigma_us_counterfactuals.csv", df_us_filtered)

## =========================================================================
##  EU Markov Switching Model
## =========================================================================
println("\n=== EU Markov Switching Model ===")

# Estimation on the same pre-2020 subsample as US for symmetry
if USE_MU_REGRESSOR && hasproperty(df, :mu_eu)
    if USE_MU_TREND
        mu_eu_full  = hp_trend(df.mu_eu, HP_LAMBDA_MO)
        mu_eu_est   = mu_eu_full[1:EST_END]
    else
        mu_eu_full  = df.mu_eu
        mu_eu_est   = df_est.mu_eu
    end
    if MU_SHARED_GAMMA
        X_sw_eu_est  = df_est.lag_lsigma_eu
        X_sw_eu_full = df.lag_lsigma_eu
        X_ns_eu_est  = reshape(collect(mu_eu_est),  :, 1)
        X_ns_eu_full = reshape(collect(mu_eu_full), :, 1)
    else
        X_sw_eu_est  = hcat(df_est.lag_lsigma_eu, mu_eu_est)
        X_sw_eu_full = hcat(df.lag_lsigma_eu,     mu_eu_full)
        X_ns_eu_est  = nothing; X_ns_eu_full = nothing
    end
else
    X_sw_eu_est  = df_est.lag_lsigma_eu
    X_sw_eu_full = df.lag_lsigma_eu
    X_ns_eu_est  = nothing; X_ns_eu_full = nothing
end
model_eu = robust_msmodel(df_est.lsigma_eu, 2, X_sw_eu_est,
                          label="EU (pre-2020)", max_retries=20, min_sigma=0.001,
                          exog_ns=X_ns_eu_est)

eu_fallback = model_eu === nothing

if eu_fallback
    println("WARNING: EU two-regime estimation failed. Falling back to single-regime AR(1).")
else
    try
        summary_msm(model_eu)
    catch e
        println("Warning: Could not display EU model summary (std errors singular)")
        println("  Coefficients β: ", model_eu.β)
        println("  Residual σ:     ", model_eu.σ)
    end
end

lsigma_eu = df.lsigma_eu

if eu_fallback
    X_eu      = [ones(N-1) lsigma_eu[1:end-1]]
    Y_eu      = lsigma_eu[2:end]
    ar1_coefs = X_eu \ Y_eu
    intercept_nor_eu = ar1_coefs[1]
    phi_nor_eu       = ar1_coefs[2]
    sigma_nor_eu     = std(Y_eu - X_eu * ar1_coefs)
    intercept_scr_eu = intercept_nor_eu
    phi_scr_eu       = phi_nor_eu
    sigma_scr_eu     = sigma_nor_eu
    P_eu_ordered     = [0.99 0.01; 0.01 0.99]
    prob_nor_eu      = fill(0.5, N)
    prob_scr_eu      = fill(0.5, N)

    println("  AR(1) intercept: ", intercept_nor_eu)
    println("  AR(1) phi:       ", phi_nor_eu)
    println("  Residual σ:      ", sigma_nor_eu)
else
    idx_nor_eu, idx_scr_eu = regime_indices(model_eu)
    println("EU regime mapping: normal = MarSwitching state $idx_nor_eu, scrambling = state $idx_scr_eu")

    # Full-sample filtered probabilities via Hamilton filter (frozen pre-2020 params)
    prob_nor_eu, prob_scr_eu = hamilton_filter_ar1(
        df.lsigma_eu, X_sw_eu_full, model_eu, idx_nor_eu, idx_scr_eu;
        X_ns = X_ns_eu_full
    )
    println("EU filtered probabilities (full sample, Hamilton): ",
            "mean prob_scr = $(round(mean(prob_scr_eu), digits=3))")

    intercept_nor_eu = model_eu.β[idx_nor_eu][1]
    phi_nor_eu       = model_eu.β[idx_nor_eu][2]
    sigma_nor_eu     = model_eu.σ[idx_nor_eu]
    intercept_scr_eu = model_eu.β[idx_scr_eu][1]
    phi_scr_eu       = model_eu.β[idx_scr_eu][2]
    sigma_scr_eu     = model_eu.σ[idx_scr_eu]

    P_eu_ordered = reorder_P(model_eu.P, idx_nor_eu, idx_scr_eu)
end

df_probs_eu = DataFrame(prob_nor = prob_nor_eu, prob_scr = prob_scr_eu)
CSV.write("data/MS_sigma_eu_prob.csv", df_probs_eu)

p4 = plot(prob_nor_eu, label="Normal", title="EU Regime Probabilities")
plot!(p4, prob_scr_eu, label="Scrambling")
savefig(p4, "data/MS_sigma_eu_probs.png")

# Simulate
if !eu_fallback
    simul_eu           = generate_msm(model_eu, 10000)
    sigma_eu_sim       = simul_eu[1]
    sigma_eu_state_sim = simul_eu[2]
else
    Random.seed!(42)
    sigma_eu_sim       = zeros(10000)
    sigma_eu_state_sim = ones(Int, 10000)
    sigma_eu_sim[1]    = mean(lsigma_eu)
    for t in 2:10000
        sigma_eu_sim[t] = intercept_nor_eu + phi_nor_eu * sigma_eu_sim[t-1] + sigma_nor_eu * randn()
    end
end
aux_eu = DataFrame(sigma_eu_sim = sigma_eu_sim, state_eu_sim = sigma_eu_state_sim)
CSV.write("data/MS_sigma_eu_simul.csv", aux_eu)

println("Estimated AR(1) coefficient (phi_eu, normal): ", phi_nor_eu)

# EU counterfactual
lsigma_eu_filtered = zeros(N)
for t in 1:N
    if t == 1
        lsigma_eu_filtered[t] = lsigma_eu[t]
    elseif prob_scr_eu[t-1] < cutoff && prob_scr_eu[t] >= cutoff
        lsigma_eu_filtered[t] = phi_nor_eu * lsigma_eu[t-1] + intercept_nor_eu
    elseif prob_scr_eu[t-1] >= cutoff && prob_scr_eu[t] >= cutoff
        lsigma_eu_filtered[t] = phi_nor_eu * lsigma_eu_filtered[t-1] + intercept_nor_eu
    else
        lsigma_eu_filtered[t] = lsigma_eu[t]
    end
end

sigma_eu_filtered = exp.(lsigma_eu_filtered)

p5 = plot(exp.(lsigma_eu), label="Actual", title="EU Sigma: Actual vs Counterfactual")
plot!(p5, sigma_eu_filtered, label="Counterfactual")
savefig(p5, "data/MS_sigma_eu_counterfactual.png")

df_eu_filtered = DataFrame(sigma_eu_filtered = sigma_eu_filtered)
CSV.write("data/MS_sigma_eu_counterfactuals.csv", df_eu_filtered)

## =========================================================================
##  Export estimated parameters for MATLAB (setup_markov.m)
## =========================================================================
println("\n=== Exporting Parameters ===")
println("US Transition matrix (row-stochastic, rows = (nor, scr)):")
println(P_us_ordered)
println("Row sums: ", sum(P_us_ordered, dims=2))

# Convention: r1 = normal, r2 = scrambling — applies to AR(1) params AND P matrix.
# When USE_MU_REGRESSOR=true, also export gamma_mu_us_r* (the AR coefficient on μ)
# and alpha_mu_us_r* = γ / (1−ρ) (the long-run μ-loading).
params_rows = [
    ("mu_sigma_us_r1",    mu_nor_us),
    ("rho_sigma_us_r1",   phi_nor_us),
    ("Sigma_sigma_us_r1", sigma_nor_us),
    ("mu_sigma_us_r2",    mu_scr_us),
    ("rho_sigma_us_r2",   phi_scr_us),
    ("Sigma_sigma_us_r2", sigma_scr_us),
    ("P11", P_us_ordered[1,1]),
    ("P12", P_us_ordered[1,2]),
    ("P21", P_us_ordered[2,1]),
    ("P22", P_us_ordered[2,2]),
]
if USE_MU_REGRESSOR && length(model_us.β[idx_nor_us]) >= 3
    gamma_nor = model_us.β[idx_nor_us][end]
    gamma_scr = model_us.β[idx_scr_us][end]
    alpha_nor = gamma_nor / (1 - phi_nor_us)
    alpha_scr = gamma_scr / (1 - phi_scr_us)
    if MU_SHARED_GAMMA
        push!(params_rows, ("gamma_mu_us_shared", gamma_nor))   # nor == scr by constraint
        push!(params_rows, ("alpha_mu_us_r1",     alpha_nor))
        push!(params_rows, ("alpha_mu_us_r2",     alpha_scr))
        println("\nμ-regressor (shared γ; α = γ/(1-ρ_k)):")
        println("  γ (shared)         = $(round(gamma_nor, digits=4))")
        println("  α (R1, normal)     = $(round(alpha_nor, digits=4))")
        println("  α (R2, scrambling) = $(round(alpha_scr, digits=4))")
    else
        push!(params_rows, ("gamma_mu_us_r1", gamma_nor))
        push!(params_rows, ("gamma_mu_us_r2", gamma_scr))
        push!(params_rows, ("alpha_mu_us_r1", alpha_nor))
        push!(params_rows, ("alpha_mu_us_r2", alpha_scr))
        println("\nμ-regressor (per-regime γ; α = γ/(1-ρ_k)):")
        println("  R1 (normal):       γ = $(round(gamma_nor, digits=4))  α = $(round(alpha_nor, digits=4))")
        println("  R2 (scrambling):   γ = $(round(gamma_scr, digits=4))  α = $(round(alpha_scr, digits=4))")
    end
end

params_df = DataFrame(
    parameter = [r[1] for r in params_rows],
    value     = [r[2] for r in params_rows],
)

CSV.write("data/MS_params.csv", params_df)

println("\nExported parameters:")
println(params_df)

println("\n=== Markov Estimation Complete ===")
println("Output files saved to data/ folder")
