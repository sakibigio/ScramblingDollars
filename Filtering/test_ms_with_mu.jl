# Test: does MarSwitching.jl accept multiple switching regressors?
# Model: log σ_t = c_k + ρ_k · log σ_{t-1} + γ_k · μ_t + ε_{k,t}
# Equivalent to user's spec with α_k = γ_k / (1 - ρ_k):
#   log σ_t = c_k + α_k·(1-ρ_k)·μ_t + ρ_k · log σ_{t-1} + ε_{k,t}

using MarSwitching
using CSV
using DataFrames
using Random
using Statistics
using LinearAlgebra

# Load RW_shock.csv (current σ_us series) and also need μ_us from data
df = CSV.read("RW_shock.csv", DataFrame, missingstring = "NA")
df.lag_sigma_us = [missing; df.sigma_us[1:end-1]]
dropmissing!(df, :lag_sigma_us)
df.lsigma_us     = log.(df.sigma_us)
df.lag_lsigma_us = log.(df.lag_sigma_us)

# Load μ_us from LFX_data.mat — but we can't load .mat from julia easily.
# Workaround: write μ_us from MATLAB into a temporary CSV first.
# Look for it on disk; if not present, instruct user.
mu_file = "_mu_us_for_test.csv"
if !isfile(mu_file)
    println("ERROR: $mu_file not found.")
    println("  From MATLAB: writematrix(mu_us, '_mu_us_for_test.csv');")
    exit(1)
end

mu_us = vec(Matrix(CSV.read(mu_file, DataFrame; header=false)))

# RW_shock spans rows 1..295 in the underlying data; after dropmissing, df has 294 rows
# corresponding to original rows 2..295. So we use mu_us[2:end].
println("Sizes — sigma_us(after drop): $(nrow(df)), mu_us(raw): $(length(mu_us))")
if length(mu_us) == nrow(df) + 1
    mu_us_aligned = mu_us[2:end]
elseif length(mu_us) == nrow(df)
    mu_us_aligned = mu_us
else
    error("mu_us length $(length(mu_us)) doesn't match df length $(nrow(df))")
end

println("\n=== Test 1: existing model (lag only) ===")
Random.seed!(0)
m1 = MSModel(df.lsigma_us, 2, exog_switching_vars = df.lag_lsigma_us)
println("β[1] = $(round.(m1.β[1], digits=4))   # state 1: [intercept, ρ]")
println("β[2] = $(round.(m1.β[2], digits=4))   # state 2: [intercept, ρ]")
println("σ    = $(round.(m1.σ, digits=4))")
println("ll = $(round(m1.Likelihood, digits=2))")

println("\n=== Test 2: with μ as additional switching regressor ===")
# Pack two regressors into a matrix: column 1 = lag log σ, column 2 = μ
X_sw = hcat(df.lag_lsigma_us, mu_us_aligned)
println("exog_switching_vars shape: $(size(X_sw))")
Random.seed!(0)
try
    m2 = MSModel(df.lsigma_us, 2, exog_switching_vars = X_sw)
    println("β[1] = $(round.(m2.β[1], digits=4))   # state 1: [intercept, ρ, γ_μ]")
    println("β[2] = $(round.(m2.β[2], digits=4))   # state 2: [intercept, ρ, γ_μ]")
    println("σ    = $(round.(m2.σ, digits=4))")
    println("ll = $(round(m2.Likelihood, digits=2))")

    # Translate to user's parameterization: α_k = γ_k / (1 - ρ_k)
    for k in 1:2
        intercept = m2.β[k][1]
        rho       = m2.β[k][2]
        gamma     = m2.β[k][3]
        alpha     = gamma / (1 - rho)
        mu_lr     = intercept / (1 - rho)
        println("State $k: intercept=$(round(intercept,digits=4)), ρ=$(round(rho,digits=4)), γ_μ=$(round(gamma,digits=4))")
        println("         => α (long-run μ-loading) = $(round(alpha,digits=4))")
        println("         => long-run mean log σ at μ=mean(μ) = $(round(mu_lr + alpha*mean(mu_us_aligned), digits=4))")
    end

    # Likelihood comparison
    println("\nLog-likelihood: lag-only = $(round(m1.Likelihood,digits=2)),  with μ = $(round(m2.Likelihood,digits=2))")
    println("Δll = $(round(m2.Likelihood - m1.Likelihood, digits=2))   (one extra parameter per state = 2 extra)")
catch e
    println("FAILED: $e")
end
