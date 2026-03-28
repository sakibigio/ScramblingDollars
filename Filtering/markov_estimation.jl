# Markov Switching Estimation for Sigma Shocks
# Part of: Scrambling for Dollars pipeline
# 
# Reads: RW_shock.csv (from main_filter.m)
# Writes: 
#   data/MS_sigma_us_prob.csv          - filtered regime probabilities
#   data/MS_sigma_us_counterfactuals.csv - counterfactual sigma path
#   data/MS_sigma_us_simul.csv         - simulated paths
#   data/MS_sigma_eu_*.csv             - same for EU
#   data/MS_params.csv                 - estimated parameters for setup_markov.m

using MarSwitching
using CSV
using DataFrames
using Plots
using Random
using Distributions
using StatsBase
using Statistics

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
    # Check if either regime σ is near zero
    if minimum(model.σ) < min_sigma
        return true, "min σ = $(minimum(model.σ)) < $min_sigma"
    end
    # Check if filtered probs put all mass on one state
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

function robust_msmodel(y, n_states, exog_sv; label="", max_retries=20, min_sigma=0.01)
    best_model = nothing
    best_ll = -Inf

    for i in 0:max_retries
        if i > 0
            Random.seed!(i)
        end
        try
            m = MSModel(y, n_states, exog_switching_vars = exog_sv)
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

    # Report final solution
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
##  US Markov Switching Model
## =========================================================================
println("\n=== US Markov Switching Model ===")

model_us = robust_msmodel(df.lsigma_us, 2, df.lag_lsigma_us, label="US")
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

# Filtered probabilities
probs_us = filtered_probs(model_us)
p2 = plot(probs_us[:,1], label="State 1", title="US Regime Probabilities")
plot!(p2, probs_us[:,2], label="State 2")
savefig(p2, "data/MS_sigma_us_probs.png")

# Save probabilities
df_probs_us = DataFrame(prob_state1 = probs_us[:,1], prob_state2 = probs_us[:,2])
CSV.write("data/MS_sigma_us_prob.csv", df_probs_us)

# Simulate
simul_us = generate_msm(model_us, 10000)
sigma_us_sim = simul_us[1]
sigma_us_state_sim = simul_us[2]
aux_us = DataFrame(sigma_us_sim = sigma_us_sim, state_us_sim = sigma_us_state_sim)
CSV.write("data/MS_sigma_us_simul.csv", aux_us)

# Counterfactuals - identify low-volatility regime
if mean(probs_us[:,1]) > mean(probs_us[:,2])
    intercept_us = model_us.β[1][1]
    phi_us = model_us.β[1][2]
    sigma_resid_us_r1 = model_us.σ[1]
    sigma_resid_us_r2 = model_us.σ[2]
    intercept_us_r2 = model_us.β[2][1]
    phi_us_r2 = model_us.β[2][2]
    state2_prob_us = probs_us[:,2]
else
    intercept_us = model_us.β[2][1]
    phi_us = model_us.β[2][2]
    sigma_resid_us_r1 = model_us.σ[2]
    sigma_resid_us_r2 = model_us.σ[1]
    intercept_us_r2 = model_us.β[1][1]
    phi_us_r2 = model_us.β[1][2]
    state2_prob_us = probs_us[:,1]
end

println("Estimated AR(1) coefficient (phi_us): ", phi_us)

# ── Export US estimated parameters to CSV ──────────────────────────────────
# Regime labels: nor = normal (high-prob state), scr = scrambling (low-prob state)
# After regime identification above:
#   nor: intercept_us, phi_us, sigma_resid_us_r1
#   scr: intercept_us_r2, phi_us_r2, sigma_resid_us_r2

mu_us_nor = intercept_us / (1 - phi_us)
mu_us_scr = intercept_us_r2 / (1 - phi_us_r2)
sigma_ss_nor = exp(mu_us_nor)
sigma_ss_scr = exp(mu_us_scr)

# Transition probabilities (leaving probabilities)
# model_us.P is the transition matrix from MarSwitching
P_us_tmp = model_us.P
if mean(probs_us[:,1]) > mean(probs_us[:,2])
    # State 1 = normal, State 2 = scrambling
    trans_nor = 1 - P_us_tmp[1,1]  # prob of leaving normal
    trans_scr = 1 - P_us_tmp[2,2]  # prob of leaving scrambling
else
    trans_nor = 1 - P_us_tmp[2,2]
    trans_scr = 1 - P_us_tmp[1,1]
end

# Standard errors via Hessian — use a function to avoid try-block scoping issues
function compute_msm_se(model, probs, ss_nor, ss_scr, phi_nor, phi_scr)
    nans = (NaN, NaN, NaN, NaN, NaN, NaN)
    try
        raw_se = MarSwitching.get_std_errors(model)
        V_σ, V_β = MarSwitching.vec2param_switch(raw_se, model.k,
                        model.n_β, model.n_β_ns, model.switching_var)
        if mean(probs[:,1]) > mean(probs[:,2])
            se_Σ_nor = V_σ[1]; se_Σ_scr = V_σ[2]
            se_int_nor = V_β[1][1]; se_int_scr = V_β[2][1]
            se_ρ_nor = V_β[1][2]; se_ρ_scr = V_β[2][2]
        else
            se_Σ_nor = V_σ[2]; se_Σ_scr = V_σ[1]
            se_int_nor = V_β[2][1]; se_int_scr = V_β[1][1]
            se_ρ_nor = V_β[2][2]; se_ρ_scr = V_β[1][2]
        end
        # Delta method: SE(exp(mu)) ≈ exp(mu) * SE(intercept) / |1 - phi|
        se_ss_nor = ss_nor * se_int_nor / abs(1 - phi_nor)
        se_ss_scr = ss_scr * se_int_scr / abs(1 - phi_scr)
        println("  Standard errors computed successfully")
        return (se_ss_nor, se_ss_scr, se_ρ_nor, se_ρ_scr, se_Σ_nor, se_Σ_scr)
    catch e
        println("  Warning: could not compute standard errors — $e")
        return nans
    end
end

(se_sigma_ss_nor, se_sigma_ss_scr, se_rho_nor, se_rho_scr, se_Sigma_nor, se_Sigma_scr) =
    compute_msm_se(model_us, probs_us, sigma_ss_nor, sigma_ss_scr, phi_us, phi_us_r2)

us_params_df = DataFrame(
    param = [
        "sigma_ss_scr",    "sigma_ss_nor",
        "sigma_ss_scr_se", "sigma_ss_nor_se",
        "rho_scr",         "rho_nor",
        "rho_scr_se",      "rho_nor_se",
        "Sigma_scr",       "Sigma_nor",
        "Sigma_scr_se",    "Sigma_nor_se",
        "trans_scr",       "trans_nor"
    ],
    value = [
        sigma_ss_scr,      sigma_ss_nor,
        se_sigma_ss_scr,   se_sigma_ss_nor,
        phi_us_r2,         phi_us,
        se_rho_scr,        se_rho_nor,
        sigma_resid_us_r2, sigma_resid_us_r1,
        se_Sigma_scr,      se_Sigma_nor,
        trans_scr,          trans_nor
    ]
)

CSV.write("data/MS_sigma_us_params.csv", us_params_df)
println("\nUS params written to data/MS_sigma_us_params.csv:")
println(us_params_df)

# Generate counterfactual (replace high-vol regimes with AR(1) from low-vol)
lsigma_us = df.lsigma_us
N = length(lsigma_us)
lsigma_us_filtered = zeros(N)
cutoff = 0.25

for t in 1:N
    if t == 1
        lsigma_us_filtered[t] = lsigma_us[t]
    elseif state2_prob_us[t-1] < cutoff && state2_prob_us[t] >= cutoff
        # Transition into high-vol state
        lsigma_us_filtered[t] = phi_us * lsigma_us[t-1] + intercept_us
    elseif state2_prob_us[t-1] >= cutoff && state2_prob_us[t] >= cutoff
        # Stay in high-vol state
        lsigma_us_filtered[t] = phi_us * lsigma_us_filtered[t-1] + intercept_us
    else
        # Low-vol state - use actual data
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

model_eu = robust_msmodel(df.lsigma_eu, 2, df.lag_lsigma_eu,
                          label="EU", max_retries=20, min_sigma=0.001)

eu_fallback = model_eu === nothing

if eu_fallback
    # Single-regime AR(1) fallback — EU series too smooth for two regimes
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
    # Single-regime AR(1) fallback
    X_eu = [ones(N-1) lsigma_eu[1:end-1]]
    Y_eu = lsigma_eu[2:end]
    ar1_coefs = X_eu \ Y_eu
    intercept_eu = ar1_coefs[1]
    phi_eu = ar1_coefs[2]
    sigma_resid_eu_r1 = std(Y_eu - X_eu * ar1_coefs)
    sigma_resid_eu_r2 = sigma_resid_eu_r1
    intercept_eu_r2 = intercept_eu
    phi_eu_r2 = phi_eu
    P_eu_mat = [0.99 0.01; 0.01 0.99]
    probs_eu = hcat(fill(0.5, N), fill(0.5, N))
    state2_prob_eu = probs_eu[:,2]

    println("  AR(1) intercept: ", intercept_eu)
    println("  AR(1) phi:       ", phi_eu)
    println("  Residual σ:      ", sigma_resid_eu_r1)
else
    # Filtered probabilities
    probs_eu = filtered_probs(model_eu)

    # Counterfactuals — identify low-volatility regime
    if mean(probs_eu[:,1]) > mean(probs_eu[:,2])
        intercept_eu = model_eu.β[1][1]
        phi_eu = model_eu.β[1][2]
        sigma_resid_eu_r1 = model_eu.σ[1]
        sigma_resid_eu_r2 = model_eu.σ[2]
        intercept_eu_r2 = model_eu.β[2][1]
        phi_eu_r2 = model_eu.β[2][2]
        state2_prob_eu = probs_eu[:,2]
    else
        intercept_eu = model_eu.β[2][1]
        phi_eu = model_eu.β[2][2]
        sigma_resid_eu_r1 = model_eu.σ[2]
        sigma_resid_eu_r2 = model_eu.σ[1]
        intercept_eu_r2 = model_eu.β[1][1]
        phi_eu_r2 = model_eu.β[1][2]
        state2_prob_eu = probs_eu[:,1]
    end
end

# Save filtered probabilities
p4 = plot(probs_eu[:,1], label="State 1", title="EU Regime Probabilities")
plot!(p4, probs_eu[:,2], label="State 2")
savefig(p4, "data/MS_sigma_eu_probs.png")

df_probs_eu = DataFrame(prob_state1 = probs_eu[:,1], prob_state2 = probs_eu[:,2])
CSV.write("data/MS_sigma_eu_prob.csv", df_probs_eu)

# Simulate
if !eu_fallback
    simul_eu = generate_msm(model_eu, 10000)
    sigma_eu_sim = simul_eu[1]
    sigma_eu_state_sim = simul_eu[2]
else
    Random.seed!(42)
    sigma_eu_sim = zeros(10000)
    sigma_eu_state_sim = ones(Int, 10000)
    sigma_eu_sim[1] = mean(lsigma_eu)
    for t in 2:10000
        sigma_eu_sim[t] = intercept_eu + phi_eu * sigma_eu_sim[t-1] + sigma_resid_eu_r1 * randn()
    end
end
aux_eu = DataFrame(sigma_eu_sim = sigma_eu_sim, state_eu_sim = sigma_eu_state_sim)
CSV.write("data/MS_sigma_eu_simul.csv", aux_eu)

println("Estimated AR(1) coefficient (phi_eu): ", phi_eu)

# Counterfactual
lsigma_eu_filtered = zeros(N)

for t in 1:N
    if t == 1
        lsigma_eu_filtered[t] = lsigma_eu[t]
    elseif state2_prob_eu[t-1] < cutoff && state2_prob_eu[t] >= cutoff
        lsigma_eu_filtered[t] = phi_eu * lsigma_eu[t-1] + intercept_eu
    elseif state2_prob_eu[t-1] >= cutoff && state2_prob_eu[t] >= cutoff
        lsigma_eu_filtered[t] = phi_eu * lsigma_eu_filtered[t-1] + intercept_eu
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

# Extract transition probabilities from model
# Try multiple approaches since MarSwitching API varies by version
function get_transition_matrix(model)
    # Try common field names
    for field in [:P, :T_mat, :trans_mat, :transition_matrix]
        if hasproperty(model, field)
            return getfield(model, field)
        end
    end
    # Try function call
    try
        return MarSwitching.transition_matrix(model)
    catch
    end
    # Last resort: print fields so user can identify
    println("Available model fields: ", fieldnames(typeof(model)))
    # Extract from raw_params — last 2 params are typically p11, p22
    rp = model.raw_params
    println("Raw params: ", rp)
    # For 2-state model with k switching vars, transition probs are last 2
    p11 = 1 / (1 + exp(-rp[end-1]))  # logit transform
    p22 = 1 / (1 + exp(-rp[end]))
    P = [p11 (1-p11); (1-p22) p22]
    println("Constructed P from raw_params (logit): ", P)
    return P
end

P_us_mat = get_transition_matrix(model_us)
println("US Transition matrix: ", P_us_mat)

if !eu_fallback
    P_eu_mat = get_transition_matrix(model_eu)
    println("EU Transition matrix: ", P_eu_mat)
end
# (If eu_fallback, P_eu_mat was already set to [0.99 0.01; 0.01 0.99] above)

# Unconditional means: mu = intercept / (1 - phi)
mu_us_r1 = intercept_us / (1 - phi_us)
mu_us_r2 = intercept_us_r2 / (1 - phi_us_r2)

println("P_us_mat BEFORE swap: ", P_us_mat)
println("Row sums: ", sum(P_us_mat, dims=2))
println("Col sums: ", sum(P_us_mat, dims=1))

if phi_us < phi_us_r2
    # Input is left-stochastic (cols sum to 1), swap + transpose to row-stochastic
    P_us_mat = [P_us_mat[2,2] 1-P_us_mat[2,2]; 1-P_us_mat[1,1] P_us_mat[1,1]]
    println("Regimes swapped to ensure R1 = normal (high persistence)")
else
    # No swap needed, but still transpose to row-stochastic
    P_us_mat = P_us_mat'
    println("No swap needed, transposed to row-stochastic")
end

println("P_us_mat AFTER: ", P_us_mat)
println("Row sums: ", sum(P_us_mat, dims=2))

# Build parameter table
# Convention: r1 = normal (high prob), r2 = volatile (low prob)
params_df = DataFrame(
    parameter = [
        "mu_sigma_us_r1", "rho_sigma_us_r1", "Sigma_sigma_us_r1",
        "mu_sigma_us_r2", "rho_sigma_us_r2", "Sigma_sigma_us_r2",
        "P11", "P12", "P21", "P22"
    ],
    value = [
        mu_us_r1, phi_us, sigma_resid_us_r1,
        mu_us_r2, phi_us_r2, sigma_resid_us_r2,
        P_us_mat[1,1], P_us_mat[1,2], P_us_mat[2,1], P_us_mat[2,2]
    ]
)

CSV.write("data/MS_params.csv", params_df)

println("\nExported parameters:")
println(params_df)
println("\nTransition matrix (US):")
println(P_us_mat)

println("\n=== Markov Estimation Complete ===")
println("Output files saved to data/ folder")
