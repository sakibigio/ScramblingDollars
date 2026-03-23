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
##  US Markov Switching Model
## =========================================================================
println("\n=== US Markov Switching Model ===")

model_us = MSModel(df.lsigma_us, 2, 
                   exog_switching_vars = df.lag_lsigma_us)

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

model_eu = MSModel(df.lsigma_eu, 2, 
                   exog_switching_vars = df.lag_lsigma_eu)

try
    summary_msm(model_eu)
catch e
    println("Warning: Could not display EU model summary (std errors singular)")
    println("  Coefficients β: ", model_eu.β)
    println("  Residual σ:     ", model_eu.σ)
end

# Filtered probabilities
probs_eu = filtered_probs(model_eu)
p4 = plot(probs_eu[:,1], label="State 1", title="EU Regime Probabilities")
plot!(p4, probs_eu[:,2], label="State 2")
savefig(p4, "data/MS_sigma_eu_probs.png")

# Save probabilities
df_probs_eu = DataFrame(prob_state1 = probs_eu[:,1], prob_state2 = probs_eu[:,2])
CSV.write("data/MS_sigma_eu_prob.csv", df_probs_eu)

# Simulate
simul_eu = generate_msm(model_eu, 10000)
sigma_eu_sim = simul_eu[1]
sigma_eu_state_sim = simul_eu[2]
aux_eu = DataFrame(sigma_eu_sim = sigma_eu_sim, state_eu_sim = sigma_eu_state_sim)
CSV.write("data/MS_sigma_eu_simul.csv", aux_eu)

# Counterfactuals
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

println("Estimated AR(1) coefficient (phi_eu): ", phi_eu)

lsigma_eu = df.lsigma_eu
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

# Transition matrix
P_us = model_us.raw_params
P_eu = model_eu.raw_params

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

P_eu_mat = get_transition_matrix(model_eu)
println("EU Transition matrix: ", P_eu_mat)

# Unconditional means: mu = intercept / (1 - phi)
mu_us_r1 = intercept_us / (1 - phi_us)
mu_us_r2 = intercept_us_r2 / (1 - phi_us_r2)

# Check if EU model is degenerate (both regimes identical)
eu_degenerate = abs(model_eu.σ[1] - model_eu.σ[2]) < 1e-6
if eu_degenerate
    println("WARNING: EU model is degenerate (regimes identical). Using single-regime AR(1).")
    # Use regime 1 params for both
    mu_eu_r1 = model_eu.β[1][1] / (1 - model_eu.β[1][2])
    mu_eu_r2 = mu_eu_r1
    phi_eu_r1_val = model_eu.β[1][2]
    phi_eu_r2_val = phi_eu_r1_val
    sigma_eu_r1_val = model_eu.σ[1]
    sigma_eu_r2_val = sigma_eu_r1_val
    # Essentially one regime
    P_eu_mat = [0.99 0.01; 0.01 0.99]
    println("  Using single-regime parameters with near-identity transition matrix")
end

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
