% filter_diagnostics.m
% Extracted from main_filter.m (lines 536–599) on 2026-02-28
% Called by: main_filter.m
%
% Contents:
%   - Steady state averages
%   - Parameter summary
%   - Sigma distribution stats
%   - Solver success rates
%   - Residual analysis (model vs data TED)
%   - Volume ratios (FF, DW) model vs data
%   - Correlation and std dev comparison (BP, TED, CIP)

%% Steady state averages
sigma_us_av = mean(sigma_us_TED_t);
sigma_eu_av = mean(sigma_eu_TED_t);
Rm_us_av = mean(Rm_us);
Rm_eu_av = mean(Rm_eu);
riskprm_av = mean(riskprm_t);
Theta_d_us_av = mean(Theta_d_us_t);
Theta_d_eu_av = mean(Theta_d_us_t);

%% Diagnostics
fprintf('\n=== Filter Diagnostics ===\n');
fprintf('lambda=%.2f, eta=%.2f, ploss=%.2f, iota_us=%.6f\n', lambda_us, eta, ploss_us, iota_us);
fprintf('sigma_us: mean=%.4f, min=%.4f, max=%.4f, p05=%.4f, p95=%.4f\n', mean(sigma_us_t), min(sigma_us_t), max(sigma_us_t), prctile(sigma_us_t,5), prctile(sigma_us_t,95));
fprintf('sigma_eu: mean=%.4f, min=%.4f, max=%.4f\n', mean(sigma_eu_t), min(sigma_eu_t), max(sigma_eu_t));
fprintf('Post/Pre(48): %.2f\n', mean(sigma_us_t(end-47:end))/mean(sigma_us_t(1:48)));
fprintf('Implied Std(omega): mean=%.1f%%, p95=%.1f%%\n', mean(sigma_us_t)*sqrt(8)*100, prctile(sigma_us_t,95)*sqrt(8)*100);

% Solver
n_fsolve = sum(sigma_us_TED_flag > 0 & sigma_us_TED_flag ~= 10);
n_fminbnd = sum(sigma_us_TED_flag == 10);
n_fail = sum(sigma_us_TED_flag <= 0);
fprintf('Solver: fsolve=%d (%.1f%%), fminbnd=%d (%.1f%%), fail=%d\n', n_fsolve, n_fsolve/T*100, n_fminbnd, n_fminbnd/T*100, n_fail);

% Residuals
model_ted = zeros(T,1); data_ted = zeros(T,1);
for tt2 = 1:T
    model_ted(tt2) = Chi_p_psi(exp(mu_us(tt2)), ploss_us, sigma_us_t(tt2), iota_us, lambda_us, eta, matching_type) * abs_scale;
    if matching_type == 0
        data_ted(tt2) = min_test_us + (TED_s_us_t(tt2) - min(TED_s_us_t)) * abs_scale;
    else
        data_ted(tt2) = TED_s_us_t(tt2) * abs_scale;
    end
end
resid = abs(model_ted - data_ted);
fprintf('Residuals(bps): max=%.2f, mean=%.2f, median=%.2f, pct<0.1=%.1f%%\n', max(resid), mean(resid), median(resid), sum(resid<0.1)/T*100);

% Volume ratios
fprintf('\n--- Volume Ratios (Model vs Data, %% of checkable deposits) ---\n');
fprintf('FF:  model=%.2f%%\n', mean(FF_us_t)*100);
fprintf('DW:  model=%.2f%%\n', mean(DW_us_t)*100);
fprintf('DW+FF=%.2f%%, DW/FF=%.2f%%\n', mean(DW_us_t+FF_us_t)*100, mean(DW_us_t./FF_us_t)*100);
if exist('DW_n', 'var') && exist('FF_n', 'var')
    dw_valid_diag = find(~isnan(DW_n));
    ff_valid_diag = find(~isnan(FF_n));
    fprintf('DW data: mean=%.3f%% (%d periods)\n', mean(DW_n(dw_valid_diag))*100, length(dw_valid_diag));
    if ~isempty(ff_valid_diag)
        fprintf('FF data: mean=%.2f%% (%d periods)\n', mean(FF_n(ff_valid_diag))*100, length(ff_valid_diag));
    end
    
    % Regime breakdown: Pre-GFC / GFC / Post-GFC
    pre  = 1:72;      % Jan 2001 – Dec 2006
    gfc  = 73:120;    % Jan 2007 – Dec 2010
    post = 121:234;   % Jan 2011 – Jun 2020
    fprintf('\n--- Volume by Regime ---\n');
    fprintf('              Pre-GFC    GFC        Post-GFC\n');
    fprintf('DW data:      %.3f%%     %.3f%%     %.3f%%\n', nanmean(DW_n(pre))*100, nanmean(DW_n(gfc))*100, nanmean(DW_n(post))*100);
    fprintf('DW model:     %.3f%%     %.3f%%     %.3f%%\n', mean(DW_us_t(pre))*100, mean(DW_us_t(gfc))*100, mean(DW_us_t(post))*100);
    fprintf('FF data:      %.2f%%     %.2f%%     %.2f%%\n', nanmean(FF_n(pre))*100, nanmean(FF_n(gfc))*100, nanmean(FF_n(post))*100);
    fprintf('FF model:     %.2f%%     %.2f%%     %.2f%%\n', mean(FF_us_t(pre))*100, mean(FF_us_t(gfc))*100, mean(FF_us_t(post))*100);
end

% Correlations and std devs
fprintf('\n--- Model vs Data ---\n');
dBP_m = diff(BP_us_t(datesperiod)); dBP_d = diff(Rb_Rm(datesperiod));
dTED_m = diff(TED_us_t(datesperiod)); dTED_d = diff(TED_s_us_t(datesperiod));
dCIP_m = diff(CIP_t(datesperiod)); dCIP_d = diff(cip(datesperiod));
fprintf('         corr(level)  corr(diff)  std_model  std_data  dstd_model  dstd_data\n');
fprintf('BP:      %.3f         %.3f        %.1f       %.1f      %.1f        %.1f\n', ...
    corr(BP_us_t(datesperiod), Rb_Rm(datesperiod)), corr(dBP_m, dBP_d), ...
    std(BP_us_t(datesperiod))*abs_scale, std(Rb_Rm(datesperiod))*abs_scale, ...
    std(dBP_m)*abs_scale, std(dBP_d)*abs_scale);
fprintf('TED:     %.3f         %.3f        %.1f       %.1f      %.1f        %.1f\n', ...
    corr(TED_us_t(datesperiod), TED_s_us_t(datesperiod)), corr(dTED_m, dTED_d), ...
    std(TED_us_t(datesperiod))*abs_scale, std(TED_s_us_t(datesperiod))*abs_scale, ...
    std(dTED_m)*abs_scale, std(dTED_d)*abs_scale);
fprintf('CIP:     %.3f         %.3f        %.1f       %.1f      %.1f        %.1f\n', ...
    corr(CIP_t(datesperiod), cip(datesperiod)), corr(dCIP_m, dCIP_d), ...
    std(CIP_t(datesperiod))*abs_scale, std(cip(datesperiod))*abs_scale, ...
    std(dCIP_m)*abs_scale, std(dCIP_d)*abs_scale);

clear n_fsolve n_fminbnd n_fail model_ted data_ted resid tt2 dBP_m dBP_d dTED_m dTED_d dCIP_m dCIP_d;
