
using MarSwitching
using CSV
using DataFrames
using Plots
using Random
using Distributions
using StatsBase
using Statistics

df = CSV.read("/Users/senaustuner/Desktop/RW_shock.csv", DataFrame, missingstring = "NA")


# Create a DataFrame with the exogenous variables
df.lag_sigma_us = [missing; df.sigma_us[1:end-1]]
df.lag_sigma_eu = [missing; df.sigma_eu[1:end-1]]
dropmissing!(df, :lag_sigma_eu) 
dropmissing!(df, :lag_sigma_us)

# Plot the data
plot([df.sigma_eu df.sigma_us])
df.lsigma_us= log.(df.sigma_us)
df.lag_lsigma_us= log.(df.lag_sigma_us)
df.lsigma_eu= log.(df.sigma_eu)
df.lag_lsigma_eu= log.(df.lag_sigma_eu)

# Run MS model for US
model = MSModel(df.lsigma_us, 2, 
                exog_switching_vars = df.lag_lsigma_us)

summary_msm(model)

plot( [filtered_probs(model)[:,1] filtered_probs(model)[:,2]])
simul_us=generate_msm(model,10000)
sigma_us_sim=simul_us[1]
sigma_us_state_sim=simul_us[2]
aux_us=DataFrame(:sigma_us_sim=>sigma_us_sim,:state_us_sim=>sigma_us_state_sim)
# Saving Outcome
CSV.write("/Users/senaustuner/Desktop/Output/MS_sigma_us_simul.csv",aux_us)

# Plot simulated data
plot(exp.(sigma_us_sim[1:100]))

# Counterfactuals
if mean(filtered_probs(model)[:,1]) > mean(filtered_probs(model)[:,2])
    intercept = model.β[1][1]
    phi_us = model.β[1][2]
    state2_prob = filtered_probs(model)[:,2]
else
    intercept = model.β[2][1]
    phi_us = model.β[2][2]
    state2_prob = filtered_probs(model)[:,1]
end

println("Estimated AR(1) coefficient (phi) : ", phi_us)

lsigma_us = df.lsigma_us
N = length(lsigma_us)
lsigma_us_filtered = zeros(N,1)

cutoff=0.25;
for t in 1:N
    if t==1
        lsigma_us_filtered[t] = lsigma_us[t]
    elseif  state2_prob[t-1] < cutoff && state2_prob[t] >= cutoff #state1-state2 case
        lsigma_us_filtered[t] = phi_us*lsigma_us[t-1] + intercept
    elseif state2_prob[t-1] >= cutoff && state2_prob[t] >= cutoff #state2-state2 case
        lsigma_us_filtered[t] = phi_us*lsigma_us_filtered[t-1] + intercept
    elseif (state2_prob[t-1] < cutoff && state2_prob[t] < cutoff) || (state2_prob[t-1] >= cutoff && state2_prob[t] < cutoff)
       lsigma_us_filtered[t]=lsigma_us[t]
    end
end

lsigma_us_filtered = vec(lsigma_us_filtered)
sigma_us = exp.(lsigma_us)
sigma_us_filtered = exp.(lsigma_us_filtered)

plot(sigma_us)
plot!(sigma_us_filtered)

df_us_filtered = DataFrame(sigma_us_filtered=sigma_us_filtered)
CSV.write("/Users/senaustuner/Desktop/Output/MS_sigma_us_counterfactuals.csv",df_us_filtered)

# Run MS model for EU
model = MSModel(df.lsigma_eu, 2, 
                exog_switching_vars = df.lag_lsigma_eu)

summary_msm(model)

plot( [filtered_probs(model)[:,1] filtered_probs(model)[:,2]])
simul_eu=generate_msm(model,10000)
sigma_eu_sim=simul_eu[1]
sigma_eu_state_sim=simul_eu[2]
aux_eu=DataFrame(:sigma_eu_sim=>sigma_eu_sim,:state_eu_sim=>sigma_eu_state_sim)
# Saving Outcome
CSV.write("/Users/senaustuner/Desktop/Output/MS_sigma_eu_simul.csv",aux_eu)

# Plot simulated data
plot(exp.(sigma_eu_sim[1:100]))

# Counterfactuals
if mean(filtered_probs(model)[:,1]) > mean(filtered_probs(model)[:,2])
    intercept = model.β[1][1]
    phi_eu = model.β[1][2]
    state2_prob = filtered_probs(model)[:,2]
else
    intercept = model.β[2][1]
    phi_eu = model.β[2][2]
    state2_prob = filtered_probs(model)[:,1]
end

println("Estimated AR(1) coefficient (phi) : ", phi_eu)

lsigma_eu = df.lsigma_eu
N = length(lsigma_eu)
lsigma_eu_filtered = zeros(N,1)

cutoff=0.25;
for t in 1:N
    if t==1
        lsigma_eu_filtered[t] = lsigma_eu[t]
    elseif  state2_prob[t-1] < cutoff && state2_prob[t] >= cutoff #state1-state2 case
        lsigma_eu_filtered[t] = phi_eu*lsigma_eu[t-1] + intercept
    elseif state2_prob[t-1] >= cutoff && state2_prob[t] >= cutoff #state2-state2 case
        lsigma_eu_filtered[t] = phi_eu*lsigma_eu_filtered[t-1] + intercept
    elseif (state2_prob[t-1] < cutoff && state2_prob[t] < cutoff) || (state2_prob[t-1] >= cutoff && state2_prob[t] < cutoff)
       lsigma_eu_filtered[t]=lsigma_eu[t]
    end
end

lsigma_eu_filtered = vec(lsigma_eu_filtered)
sigma_eu = exp.(lsigma_eu)
sigma_eu_filtered = exp.(lsigma_eu_filtered)

plot(sigma_eu)
plot!(sigma_eu_filtered)

df_eu_filtered = DataFrame(sigma_eu_filtered=sigma_eu_filtered)
CSV.write("/Users/senaustuner/Desktop/Output/MS_sigma_eu_factuals.csv",df_eu_filtered)