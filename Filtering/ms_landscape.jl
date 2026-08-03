# ms_landscape.jl
#
# Map the local-optimum landscape of the 2-regime Markov-switching AR(1) MLE
# on the eta=0.5 filtered sigma series (RW_shock_eta50.csv), truncated pre-2020.
#
# Starts: 50 random seeds (k-means init), a few with random_search_em, plus
# structured x0's taken from the BASELINE (cbase) series fit — raw and with
# intercepts shifted to the eta50 series' level. Reports every distinct mode
# with its log-likelihood so "flat likelihood" can be judged directly.
#
# Standalone: reads tagged CSVs only, writes ms_landscape_results.csv.

using MarSwitching, CSV, DataFrames, Random, Statistics

const EST_END = 227

function load_series(path)
    df = CSV.read(path, DataFrame)
    ls = log.(df.sigma_us)
    y   = ls[2:end]
    lag = ls[1:end-1]
    y[1:EST_END], reshape(lag[1:EST_END], :, 1)
end

# regime summary: sort regimes so regime 1 = lower Sigma (the "normal" one)
function summarize(m)
    k = 2
    rho = [m.β[i][2] for i in 1:k]
    ic  = [m.β[i][1] for i in 1:k]
    ord = sortperm(m.σ)                       # low variance first = normal
    ic, rho, sg = ic[ord], rho[ord], m.σ[ord]
    Pd = [m.P[ord[1],ord[1]], m.P[ord[2],ord[2]]]
    lr = [abs(1 - r) < 1e-9 ? NaN : i / (1 - r) for (i, r) in zip(ic, rho)]
    (ll = m.Likelihood, mu_nor = lr[1], mu_scr = lr[2], sep = lr[2] - lr[1],
     rho_nor = rho[1], rho_scr = rho[2], sig_nor = sg[1], sig_scr = sg[2],
     P11 = Pd[1], P22 = Pd[2])
end

y50, x50   = load_series("RW_shock_eta50.csv")
ybase, xbase = load_series("RW_shock_cbase.csv")

println("eta50 series: n=$(length(y50)) mean=$(round(mean(y50), digits=3))")
println("cbase series: n=$(length(ybase)) mean=$(round(mean(ybase), digits=3))")

# --- reference fit on the baseline series (best of 5 seeds) ---
best_base = nothing
for s in 0:4
    Random.seed!(s)
    try
        m = MSModel(ybase, 2, exog_switching_vars = xbase)
        if !isnan(m.Likelihood) && (best_base === nothing || m.Likelihood > best_base.Likelihood)
            global best_base = m
        end
    catch; end
end
sb = summarize(best_base)
println("\ncbase reference: ll=$(round(sb.ll, digits=2)) mu_nor=$(round(sb.mu_nor, digits=3)) mu_scr=$(round(sb.mu_scr, digits=3)) sep=$(round(sb.sep, digits=3))")

# --- structured x0's from the baseline fit ---
# raw_params layout (k=2, switching var, switching intercept, 1 switching beta):
#   [sig1, sig2, ic1, ic2, rho1, rho2, p_raw...]
raw = copy(best_base.raw_params)
shift = mean(y50) - mean(ybase)
raw_shifted = copy(raw)
for i in 1:2
    rho_i = raw[4 + i]
    raw_shifted[2 + i] = raw[2 + i] + shift * (1 - rho_i)
end

starts = Vector{Tuple{String, Union{Nothing, Vector{Float64}}, Int}}()
push!(starts, ("x0=cbase_raw", raw, 0))
push!(starts, ("x0=cbase_shifted", raw_shifted, 0))
for s in 0:49
    push!(starts, ("seed$s", nothing, s))
end
for s in 0:4
    push!(starts, ("seed$(s)_rsem10", nothing, s))   # with random_search_em
end

rows = DataFrame(source=String[], ll=Float64[], mu_nor=Float64[], mu_scr=Float64[],
                 sep=Float64[], rho_nor=Float64[], rho_scr=Float64[],
                 sig_nor=Float64[], sig_scr=Float64[], P11=Float64[], P22=Float64[])
for (name, x0, s) in starts
    Random.seed!(s)
    try
        m = if x0 !== nothing
            MSModel(y50, 2, exog_switching_vars = x50, x0 = x0)
        elseif endswith(name, "_rsem10")
            MSModel(y50, 2, exog_switching_vars = x50, random_search_em = 10)
        else
            MSModel(y50, 2, exog_switching_vars = x50)
        end
        isnan(m.Likelihood) && continue
        r = summarize(m)
        push!(rows, (name, r.ll, r.mu_nor, r.mu_scr, r.sep, r.rho_nor, r.rho_scr,
                     r.sig_nor, r.sig_scr, r.P11, r.P22))
    catch e
        println("  $name failed: $e")
    end
end

# distinct modes: round ll to 2dp and sep to 3dp
rows.mode = string.(round.(rows.ll, digits=2), "_", round.(rows.sep, digits=2))
gd = combine(groupby(rows, :mode),
    :ll => maximum => :ll, :mu_nor => first => :mu_nor, :mu_scr => first => :mu_scr,
    :sep => first => :sep, :rho_nor => first => :rho_nor, :rho_scr => first => :rho_scr,
    :sig_nor => first => :sig_nor, :sig_scr => first => :sig_scr,
    :P11 => first => :P11, :P22 => first => :P22, nrow => :hits,
    :source => (x -> join(unique(first.(split.(x, "_"), 1))[1:min(3,end)], ",")) => :sources)
sort!(gd, :ll, rev = true)
gd.dll = gd.ll .- maximum(gd.ll)

println("\n===== DISTINCT MODES on eta50 series (truncated, sorted by loglik) =====")
println("   ll      dll    mu_nor  mu_scr   sep    rho_n  rho_s  sig_n  sig_s   P11    P22   hits")
for r in eachrow(gd)
    println(join([lpad(round(r.ll, digits=2), 8), lpad(round(r.dll, digits=2), 7),
        lpad(round(r.mu_nor, digits=3), 8), lpad(round(r.mu_scr, digits=3), 8),
        lpad(round(r.sep, digits=3), 7), lpad(round(r.rho_nor, digits=3), 7),
        lpad(round(r.rho_scr, digits=3), 7), lpad(round(r.sig_nor, digits=3), 7),
        lpad(round(r.sig_scr, digits=3), 7), lpad(round(r.P11, digits=3), 7),
        lpad(round(r.P22, digits=3), 7), lpad(r.hits, 5)], " "), "  [", r.sources, "]")
end
CSV.write("ms_landscape_results.csv", rows)
println("\nwrote ms_landscape_results.csv ($(nrow(rows)) fits, $(nrow(gd)) distinct modes)")
