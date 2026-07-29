%% build_iota_corridor_monthly.m
% Construct a MONTHLY time-varying iota "corridor" spread from the weekly
% rates file, aligned to the model's monthly sample grid.
%
%   iota counterpart = DW - RRP  =  DWR (col 7) - RRPONTSYAWARD (col 21)
%   (the ON reverse-repo award rate as the safe-rate counterpart).
%   NOTE: RRP = 0 before 25-Sep-2013 (facility did not exist), so pre-2013
%   this equals the raw discount-window rate level (user's literal choice).
%
% Weekly -> monthly by calendar-month averaging.  Aligned onto the model
% grid datenum(2001, 1:T_model, 1).  Months with no weekly coverage are
% left NaN here and padded by the caller (2001-2002 lead, 2024-2025 tail).
%
% Output (saved to data/iota_corridor_monthly.mat):
%   iota_dates_m    T_model x 1 datenum grid (Jan-2001 .. )
%   iota_sprd_pct   T_model x 1 monthly DW-Tbill spread in PERCENT (annualized)
%   iota_sprd_dec   T_model x 1 same, as a decimal fraction (pct/100)
%   cover_mask      T_model x 1 logical: true where weekly data existed
%   (also DW-IOR and FF-IOR monthly, for reference/robustness)

F = 'data/data_Saki_july_2026_rates_v2.xlsx';

%% --- Read weekly sheet ---
R = readtable(F,'Sheet','Sheet2','Range','A1:AA100000','VariableNamingRule','preserve');
wk_dates = R{:,1};
if ~isdatetime(wk_dates), wk_dates = datetime(wk_dates,'ConvertFrom','excel'); end
good = ~isnat(wk_dates);
R = R(good,:); wk_dates = wk_dates(good);

DWR      = R{:,7};    % discount window (primary credit) rate, % annual
IORB     = R{:,9};    % interest on reserve balances, % annual (0 before Oct-2008)
RRP      = R{:,21};   % ON reverse-repo award rate, % annual (0 before Sep-2013)
Tbill_1m = R{:,22};   % 1-month T-bill (Tbond_1_month), % annual
DW_IOR   = R{:,16};   % DW - IOR (reference)
FF_IOR   = R{:,17};   % FF - IOR (reference)

%% --- SPLICED base rate: a floor that is DEFINED throughout the sample -------
%   base_t =  T-bill        for  t <  Oct-2008   (RRP & IOR = 0 then)
%             IORB          for  Oct-2008 <= t < Sep-2013
%             RRP           for  t >= Sep-2013
% This avoids the RRP=0 pre-2013 degeneracy (which made DW-RRP and the TED
% adjustment collapse to the risk-free RATE LEVEL and inflate pre-2008 TED/BP).
% base_mode: 'spliced' (Tbill->IOR->RRP) or 'tbill' (T-bill throughout).
if ~exist('base_mode','var'), base_mode = 'spliced'; end
base = Tbill_1m;                                              % pre-2008 / tbill-mode default
if strcmp(base_mode,'spliced')
    mid  = wk_dates >= datetime(2008,10,1) & wk_dates < datetime(2013,9,25);
    lat  = wk_dates >= datetime(2013,9,25);
    base(mid) = IORB(mid);
    base(lat) = RRP(lat);
else
    mid = false(size(wk_dates)); lat = false(size(wk_dates));  % all T-bill
end
fprintf('Base mode: %s\n', base_mode);

DW_Tbill = DWR - base;      % iota counterpart: DW - spliced base (percent)
% TED base-rate adjustment to ADD to the data TED (and BP): (Tbill - base).
%   pre-2008 base=Tbill  -> adj = 0 (TED stays the standard LIBOR-Tbill spread)
%   2008-13  base=IOR    -> adj = Tbill - IOR
%   2013+    base=RRP    -> adj = Tbill - RRP
TedAdj_wk = Tbill_1m - base;

fprintf('Weekly rows: %d  (%s .. %s)\n', numel(wk_dates), ...
    datestr(min(wk_dates)), datestr(max(wk_dates)));
fprintf('Spliced base: Tbill(pre-Oct08)->IOR(08-13)->RRP(13+).  n_mid=%d n_late=%d\n', sum(mid), sum(lat));
fprintf('DW-base weekly: mean=%.3f min=%.3f max=%.3f  (NaN=%d)\n', ...
    mean(DW_Tbill,'omitnan'), min(DW_Tbill), max(DW_Tbill), sum(isnan(DW_Tbill)));
fprintf('TED adj weekly (Tbill-base): mean=%.3f min=%.3f max=%.3f\n', ...
    mean(TedAdj_wk,'omitnan'), min(TedAdj_wk), max(TedAdj_wk));

%% --- Weekly -> monthly (calendar-month average) ---
wk_ym = 100*year(wk_dates) + month(wk_dates);
[uym, ~, g] = unique(wk_ym);
mon_DWTb  = accumarray(g, DW_Tbill, [], @(v) mean(v,'omitnan'));
mon_DWIOR = accumarray(g, DW_IOR,   [], @(v) mean(v,'omitnan'));
mon_FFIOR = accumarray(g, FF_IOR,   [], @(v) mean(v,'omitnan'));
% TED base-rate adjustment: (Tbill - spliced base), added to data TED & BP.
mon_TedAdj = accumarray(g, TedAdj_wk, [], @(v) mean(v,'omitnan'));
mon_year  = floor(uym/100);
mon_mon   = mod(uym,100);
mon_datenum = datenum(mon_year, mon_mon, 1);

%% --- Model monthly grid: Jan-2001 .. (match main_filter.m) ---
% length set to cover through Jul-2025 (T=295); caller can pass T_model.
if ~exist('T_model','var'), T_model = 295; end
iota_dates_m = datenum(2001, 1:T_model, 1)';

iota_sprd_pct = NaN(T_model,1);
DWIOR_m       = NaN(T_model,1);
FFIOR_m       = NaN(T_model,1);
ted_adj_pct   = NaN(T_model,1);   % (Tbond - RRP) to add to data TED
cover_mask    = false(T_model,1);
for k = 1:T_model
    j = find(mon_datenum == iota_dates_m(k), 1);
    if ~isempty(j)
        iota_sprd_pct(k) = mon_DWTb(j);
        DWIOR_m(k)       = mon_DWIOR(j);
        FFIOR_m(k)       = mon_FFIOR(j);
        ted_adj_pct(k)   = mon_TedAdj(j);
        cover_mask(k)    = ~isnan(mon_DWTb(j));
    end
end
iota_sprd_dec = iota_sprd_pct / 100;
% Convert TED adjustment to the model's per-period decimal (annual %/100/12).
ted_adj_dec   = ted_adj_pct / 100 / 12;

%% --- Report coverage & pre/post stats ---
first_cov = find(cover_mask,1); last_cov = find(cover_mask,1,'last');
fprintf('\nModel grid: %s .. %s (T=%d)\n', ...
    datestr(iota_dates_m(1),'mmm-yyyy'), datestr(iota_dates_m(end),'mmm-yyyy'), T_model);
fprintf('Weekly-data coverage on grid: %s .. %s (%d of %d months)\n', ...
    datestr(iota_dates_m(first_cov),'mmm-yyyy'), datestr(iota_dates_m(last_cov),'mmm-yyyy'), ...
    sum(cover_mask), T_model);
fprintf('Uncovered (need padding): lead=%d months (pre-%s), tail=%d months (post-%s)\n', ...
    first_cov-1, datestr(iota_dates_m(first_cov),'mmm-yyyy'), ...
    T_model-last_cov, datestr(iota_dates_m(last_cov),'mmm-yyyy'));

t_break = 94;  % Oct-2008
pre  = (1:T_model)' <= t_break & cover_mask;
post = (1:T_model)' >  t_break & cover_mask;
fprintf('\nMonthly DW-base (percent):  PRE mean=%.3f [%.3f, %.3f]   POST mean=%.3f [%.3f, %.3f]\n', ...
    mean(iota_sprd_pct(pre),'omitnan'), min(iota_sprd_pct(pre)), max(iota_sprd_pct(pre)), ...
    mean(iota_sprd_pct(post),'omitnan'), min(iota_sprd_pct(post)), max(iota_sprd_pct(post)));
fprintf('  (reference) DW-IOR pre=%.3f post=%.3f ; FF-IOR pre=%.3f post=%.3f\n', ...
    mean(DWIOR_m(pre),'omitnan'), mean(DWIOR_m(post),'omitnan'), ...
    mean(FFIOR_m(pre),'omitnan'), mean(FFIOR_m(post),'omitnan'));

%% --- Save ---
fprintf('TED adj (Tbond-RRP, pct): PRE mean=%.3f  POST mean=%.3f  (added to data TED)\n', ...
    mean(ted_adj_pct(pre),'omitnan'), mean(ted_adj_pct(post),'omitnan'));

save('data/iota_corridor_monthly.mat', ...
    'iota_dates_m', 'iota_sprd_pct', 'iota_sprd_dec', 'cover_mask', ...
    'DWIOR_m', 'FFIOR_m', 'ted_adj_pct', 'ted_adj_dec', 't_break');
fprintf('\n[saved] data/iota_corridor_monthly.mat\n');
