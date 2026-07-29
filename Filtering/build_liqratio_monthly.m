%% build_liqratio_monthly.m
% Build a MONTHLY US liquidity ratio from data_Saki_july_2026.xlsx and align
% it to the model grid, to replace mu_us (= log liquidity ratio).
%
%   Source: Sheet2, col 23 "Liquidity Ratio H_18" (a LEVEL, ~0.20-0.49).
%   Model mu_us is the LOG ratio, so mu_us_new = log(monthly-mean level).
%
% Weekly -> monthly by calendar-month average of the level.  Aligned to
% datenum(2001,1:T,1); lead/tail gaps (2001-02, 2024-25) left NaN and padded
% by the caller (hold-first / hold-last).
%
% Output: data/liq_ratio_monthly.mat
%   liq_dates_m   T x 1 datenum grid
%   mu_us_new     T x 1 log liquidity ratio (NaN where uncovered)
%   liq_level_m   T x 1 monthly level (for reference)
%   cover_mask    T x 1 logical

F = 'data/data_Saki_july_2026.xlsx';
R = readtable(F,'Sheet','Sheet2','Range','A1:AC100000','VariableNamingRule','preserve');
d = R{:,1}; if ~isdatetime(d), d = datetime(d,'ConvertFrom','excel'); end
good = ~isnat(d); R=R(good,:); d=d(good);

liq_level = R{:,23};   % Liquidity Ratio H_18 (level)
fprintf('Weekly liq ratio (H_18): %d rows (%s..%s), level mean=%.3f [%.3f,%.3f]\n', ...
    numel(d), datestr(min(d)), datestr(max(d)), mean(liq_level,'omitnan'), min(liq_level), max(liq_level));

% Weekly -> monthly (average the level within each calendar month)
ym = 100*year(d)+month(d);
[uym,~,g] = unique(ym);
mon_lvl = accumarray(g, liq_level, [], @(v) mean(v,'omitnan'));
mon_dn  = datenum(floor(uym/100), mod(uym,100), 1);

if ~exist('T_model','var'), T_model = 295; end
liq_dates_m = datenum(2001,1:T_model,1)';
liq_level_m = NaN(T_model,1);
cover_mask  = false(T_model,1);
for k = 1:T_model
    j = find(mon_dn == liq_dates_m(k),1);
    if ~isempty(j)
        liq_level_m(k) = mon_lvl(j);
        cover_mask(k)  = ~isnan(mon_lvl(j));
    end
end
mu_us_new = log(liq_level_m);

fc = find(cover_mask,1); lc = find(cover_mask,1,'last');
fprintf('Coverage: %s..%s (%d/%d months); pad lead=%d, tail=%d\n', ...
    datestr(liq_dates_m(fc),'mmm-yyyy'), datestr(liq_dates_m(lc),'mmm-yyyy'), ...
    sum(cover_mask), T_model, fc-1, T_model-lc);

t_break = 94;
pre = (1:T_model)'<=t_break & cover_mask; post = (1:T_model)'>t_break & cover_mask;
fprintf('mu_us_new (log): PRE mean=%.3f [%.3f,%.3f]  POST mean=%.3f [%.3f,%.3f]\n', ...
    mean(mu_us_new(pre),'omitnan'),min(mu_us_new(pre)),max(mu_us_new(pre)), ...
    mean(mu_us_new(post),'omitnan'),min(mu_us_new(post)),max(mu_us_new(post)));

save('data/liq_ratio_monthly.mat','liq_dates_m','mu_us_new','liq_level_m','cover_mask','t_break');
fprintf('[saved] data/liq_ratio_monthly.mat\n');
