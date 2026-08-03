function obj = obj_lambda_baseline(lam, eta_fix, iot, ploss_fix, in, tg)
% OBJ_LAMBDA_BASELINE  Baseline estimation criterion restricted to lambda.
%
% Used by the robustness grids (estimate_robustness_iota.m / _ploss.m) so
% that the re-estimated lambda* minimizes the SAME criterion as the
% baseline calibration (estimate_obj_2d in estimate_params_4d.m) with
% (iota, eta, ploss) held fixed:
%
%   obj = w_bpus    * mean(r_bpus.^2) / v_bpus     (US bond premium, demeaned, bps)
%       + w_bpeu    * mean(r_bpeu.^2) / v_bpeu     (EU bond premium, demeaned, bps)
%       + w_meanFF  * (log mean FF_m  - log mean FF_d)^2
%       + w_meanDWS * (log mean DWS_m - log mean DW_d)^2   (stock vs stock)
%       + w_corr    * corr(log sigma_us, mu_us)^2
%
% Weights in tg must match the baseline run (w_cip = 0 there, so CIP is
% omitted).  These are exactly the objects plotted in the Appendix I
% robustness figures.
%
% tg fields: idx_bpus, idx_bpeu, bpus_d_dm, bpeu_d_dm, v_bpus, v_bpeu,
%            idx_ff, idx_dw, logmean_FF_d, logmean_DW_d,
%            w_bpus, w_bpeu, w_meanFF, w_meanDWS, w_corr

    % Bounds mirror b_lambda in estimate_params_4d.m
    pen = 0;
    if lam <= 0.5,  pen = pen + 1e6 * (1 + abs(0.5 - lam)); end
    if lam >= 12.0, pen = pen + 1e6 * (1 + abs(lam - 12.0)); end
    if pen > 0, obj = pen; return; end

    try
        out = run_filter_series(lam, eta_fix, iot, ploss_fix, in);
    catch
        obj = 1e9; return;
    end

    BPus = out.BP_us_t;  BPeu = out.BP_eu_t;
    FF   = out.FF_us_t;  DW   = out.DW_us_t;  DWS = out.DWS_us_t;
    sig  = out.sigma_us_t;

    % Reject if ANY period is bad (mirrors the baseline objective's
    % full-sample feasibility check).
    if any(~isreal(BPus)) || any(~isreal(BPeu)) || ...
       any(~isfinite(BPus)) || any(~isfinite(BPeu)) || ...
       any(~isfinite(FF)) || any(~isfinite(DW)) || ...
       any(DW <= 0) || any(FF <= 0)
        obj = 1e8; return;
    end

    sc = in.abs_scale;

    % --- Price residuals (bps, demeaned over the valid data sample) ---
    bpus_m_dm = BPus(tg.idx_bpus) - mean(BPus(tg.idx_bpus));
    r_bpus    = (bpus_m_dm - tg.bpus_d_dm) * sc;
    bpeu_m_dm = BPeu(tg.idx_bpeu) - mean(BPeu(tg.idx_bpeu));
    r_bpeu    = (bpeu_m_dm - tg.bpeu_d_dm) * sc;

    % --- Volume log-mean moments (FF flow, DW stock) ---
    m_FF = log(mean(FF(tg.idx_ff))) - tg.logmean_FF_d;
    if any(DWS(tg.idx_dw) <= 0), obj = 1e8; return; end
    m_DWS = log(mean(DWS(tg.idx_dw))) - tg.logmean_DW_d;

    % --- Orthogonality penalty: corr(log sigma_us, mu_us)^2 ---
    valid_sm = isfinite(sig) & sig > 0 & isfinite(in.mu_us);
    if sum(valid_sm) > 10
        lsig = log(sig(valid_sm));
        mvec = in.mu_us(valid_sm);
        c = corr(lsig(:), mvec(:));
        if ~isfinite(c), c = 1; end
        corr_pen = c^2;
    else
        corr_pen = 1;
    end

    obj = tg.w_bpus    * mean(r_bpus.^2) / tg.v_bpus ...
        + tg.w_bpeu    * mean(r_bpeu.^2) / tg.v_bpeu ...
        + tg.w_meanFF  * m_FF^2 ...
        + tg.w_meanDWS * m_DWS^2 ...
        + tg.w_corr    * corr_pen;
end
