# ms_common_rho.jl
#
# Common-rho 2-regime Markov spec: regimes switch in INTERCEPT and VARIANCE,
# single common AR(1) coefficient (lag passed as non-switching exog_vars).
# Motivation: the switching-rho spec has no stationary optimum on the full
# sample (flat post-2020 tail absorbs one regime's rho -> 1). With a common
# rho, the tail cannot hijack persistence; test whether full-sample
# estimation becomes viable.
#
# Standalone: reads tagged RW_shock_est_*.csv only; touches no live files.

using MarSwitching, CSV, DataFrames, Random, Statistics

const TRUNC = 227

function load_series(path, nend)
    df = CSV.read(path, DataFrame)
    ls = log.(df.sigma_us)
    y   = ls[2:end]
    lag = ls[1:end-1]
    n = nend <= 0 ? length(y) : min(nend, length(y))
    y[1:n], reshape(lag[1:n], :, 1)
end

function fit_cell(name, path, nend)
    y, x = load_series(path, nend)
    best = nothing
    modes = Dict{String, Any}()
    for s in 0:24
        Random.seed!(s)
        try
            m = MSModel(y, 2, exog_vars = x)   # lag NON-switching -> common rho
            isnan(m.Likelihood) && continue
            rho = m.β[1][end]                   # common across regimes
            ic  = [m.β[i][1] for i in 1:2]
            ord = sortperm(m.σ)                 # low variance = normal
            ic  = ic[ord]; sg = m.σ[ord]
            Pd  = [m.P[ord[1],ord[1]], m.P[ord[2],ord[2]]]
            lr  = abs(1 - rho) < 1e-9 ? [NaN, NaN] : ic ./ (1 - rho)
            key = string(round(m.Likelihood, digits=1), "_", round(lr[2]-lr[1], digits=2))
            rec = (ll=m.Likelihood, rho=rho, mu_nor=lr[1], mu_scr=lr[2],
                   sep=lr[2]-lr[1], sig_nor=sg[1], sig_scr=sg[2], P11=Pd[1], P22=Pd[2])
            if !haskey(modes, key) || rec.ll > modes[key].ll
                modes[key] = rec
            end
            if best === nothing || m.Likelihood > best.ll
                best = rec
            end
        catch; end
    end
    println("\n=== $name (n=$(length(y))) ===")
    if best === nothing
        println("  all seeds failed")
        return
    end
    ms = sort(collect(values(modes)), by = r -> -r.ll)
    println("   ll      rho    mu_nor  mu_scr   sep    sig_n   sig_s   P11    P22")
    for r in ms[1:min(3, end)]
        println(join([lpad(round(r.ll, digits=1), 8), lpad(round(r.rho, digits=4), 8),
            lpad(round(r.mu_nor, digits=3), 8), lpad(round(r.mu_scr, digits=3), 8),
            lpad(round(r.sep, digits=3), 7), lpad(round(r.sig_nor, digits=3), 7),
            lpad(round(r.sig_scr, digits=3), 7), lpad(round(r.P11, digits=3), 7),
            lpad(round(r.P22, digits=3), 7)], " "))
    end
end

for (tag, path) in [("nopen_lcr", "RW_shock_est_nopen_lcr.csv"),
                    ("eta50_lcr", "RW_shock_est_eta50_lcr.csv")]
    fit_cell("$tag  TRUNCATED (pre-2020)", path, TRUNC)
    fit_cell("$tag  FULL SAMPLE",          path, 0)
end
