# Markov Switching Estimation for Sigma Shocks
# Part of: Scrambling for Dollars pipeline
# 
# Reads: RW_shock.csv (from main_filter.m)
# Writes: data/MS_sigma_us_prob.csv, data/MS_sigma_us_counterfactuals.csv, etc.

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

summary_msm(model_us)

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
    state2_prob_us = probs_us[:,2]
else
    intercept_us = model_us.β[2][1]
    phi_us = model_us.β[2][2]
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

summary_msm(model_eu)

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
    state2_prob_eu = probs_eu[:,2]
else
    intercept_eu = model_eu.β[2][1]
    phi_eu = model_eu.β[2][2]
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

println("\n=== Markov Estimation Complete ===")
println("Output files saved to data/ folder")
