%% Diagnostic: empirically verify the TED asymptote at sigma -> infty
cd('/Users/sakibigio/Dropbox/Scrabmling for Dollars/Code 2025 ScramblingDollars/Filtering');
addpath('functions');
addpath('functions/chi');
addpath('data');

freq      = 12;
abs_sc    = freq * 1e4;

% Parameters (Run 6 optimum)
iota_ss   = 0.0417;
eta       = 0.300;
varrho    = 0.0;
lambda    = 3.40;
ploss     = 0.5;
matching_type = 1;
pi_ss     = 1;

iota_loc  = iota_ss / freq / pi_ss;

% Theoretical asymptote (1-eta)*iota in different units:
asy_per_period_decimal = (1 - eta) * iota_loc;
asy_per_period_bps     = asy_per_period_decimal * 1e4;
asy_annual_decimal     = (1 - eta) * iota_ss;
asy_annual_bps         = asy_annual_decimal * 1e4;
asy_chi_times_absc     = asy_per_period_decimal * abs_sc;

fprintf('Theoretical (1-eta)*iota interpretations:\n');
fprintf('  (1-eta)*iota_loc                       = %.6e [per-period, decimal]\n', asy_per_period_decimal);
fprintf('  (1-eta)*iota_loc * 1e4                 = %.4f bps  [per-period bps]\n', asy_per_period_bps);
fprintf('  (1-eta)*iota_ss  (annual decimal)      = %.6e\n', asy_annual_decimal);
fprintf('  (1-eta)*iota_ss * 1e4 (annual bps)     = %.4f bps  ← what the constraint uses\n', asy_annual_bps);
fprintf('  (1-eta)*iota_loc * abs_sc (=12*1e4)    = %.4f bps  ← should match Chi_p_psi*abs_sc\n', asy_chi_times_absc);
fprintf('\n');

% Empirical: evaluate Chi_p_psi at very large sigma, sweep mu
mu_vec = [0.01, 0.1, 0.3, 0.5, 1.0];
fprintf('Empirical Chi_p_psi(mu, ploss, sigma=1e6, ...) * abs_sc:\n');
fprintf('   mu      Chi*abs_sc (bps)\n');
for k = 1:numel(mu_vec)
    mu  = mu_vec(k);
    val = Chi_p_psi(mu, ploss, 1e6, iota_loc, lambda, eta, matching_type, varrho);
    fprintf('  %5.3f   %10.4f\n', mu, val * abs_sc);
end
fprintf('\n');

% Sweep sigma at a typical mu
mu = 0.2;
sig_vec = [0.1, 0.5, 1, 5, 10, 15, 100, 1e4, 1e6];
fprintf('Chi_p_psi(mu=%.2f, ploss=%.2f, ..., varrho=%.2f) * abs_sc, sweeping sigma:\n', mu, ploss, varrho);
fprintf('    sigma            Chi*abs_sc (bps)\n');
for k = 1:numel(sig_vec)
    val = Chi_p_psi(mu, ploss, sig_vec(k), iota_loc, lambda, eta, matching_type, varrho);
    fprintf('  %10.4g   %12.4f\n', sig_vec(k), val * abs_sc);
end
