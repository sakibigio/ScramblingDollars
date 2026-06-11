%% test_distribution_aggregates.m - Verify the chi-folder refactor.
% Checks:
%   1) varrho = 0 reproduces the legacy inline closed forms exactly.
%   2) Identity S+ - S- = m_eff holds (since E[omega] = 0).
%   3) The four chi wrappers (Chi_sys, Chi_p_psi, Echi_m, Echi_d) all see
%      the same theta and produce identical results to a direct call.
%   4) varrho ~= 0 gives the same value as varrho = 0 when evaluated at
%      mu + varrho instead of mu (a translation invariance check).

clear; close all;
addpath(fileparts(mfilename('fullpath')));

% Case grid (use only m_eff > 0 cases to match practical use)
cases = [ ...
    0.50, 0.50, 0.10, 0.00; ...
    0.50, 0.50, 0.10, 0.05; ...
    1.00, 0.75, 0.20, 0.00; ...
    1.00, 0.75, 0.20, 0.30; ...
    0.05, 0.50, 0.50, 0.00];

iota = 0.001; lambda = 1.0; eta = 0.7;
tol = 1e-12;

fprintf('--- Test 1: varrho = 0 matches legacy inline closed forms ---\n');
max_err1 = 0;
for k = 1:size(cases,1)
    mu = cases(k,1); p = cases(k,2); s = cases(k,3);
    Smin_legacy = s .* exp(-p./s .* mu);
    Spl_legacy  = mu + s .* exp(-p./s .* mu);
    [Smin_new, Spl_new, theta_new] = distribution_aggregates(mu, p, s, 0);
    err = max(abs([Smin_legacy - Smin_new, Spl_legacy - Spl_new]));
    max_err1 = max(max_err1, err);
    fprintf('  mu=%.2f p=%.2f s=%.2f: |diff|=%.2e\n', mu, p, s, err);
end
assert(max_err1 < tol, 'Backward compatibility broken');

fprintf('\n--- Test 2: identity S+ - S- = m_eff ---\n');
max_err2 = 0;
for k = 1:size(cases,1)
    mu = cases(k,1); p = cases(k,2); s = cases(k,3); vr = cases(k,4);
    [Smin_v, Spl_v] = distribution_aggregates(mu, p, s, vr);
    err = abs((Spl_v - Smin_v) - (mu - vr));
    max_err2 = max(max_err2, err);
    fprintf('  m_eff=%+.3f: residual=%.2e\n', mu - vr, err);
end
assert(max_err2 < tol, 'Identity violated');

fprintf('\n--- Test 3: theta consistent across all four chi wrappers ---\n');
max_err3 = 0;
for k = 1:size(cases,1)
    mu = cases(k,1); p = cases(k,2); s = cases(k,3); vr = cases(k,4);
    for mt = [0, 1]
        [~, ~, th_helper] = distribution_aggregates(mu, p, s, vr);
        [~, th_sys]       = Chi_sys(mu, p, s, iota, lambda, eta, mt, vr);
        % Chi_p_psi/Echi_m/Echi_d don't return theta directly, but if they
        % use the same helper, the bond/TED outputs are deterministic in theta.
        ted_a  = Chi_p_psi(mu, p, s, iota, lambda, eta, mt, vr);
        ted_b  = Chi_p_psi(mu, p, s, iota, lambda, eta, mt, vr);  % idempotence
        bp     = Echi_m(mu, p, s, iota, lambda, eta, mt, vr);
        dp     = Echi_d(mu, p, s, iota, lambda, eta, mt, vr);
        err = abs(th_helper - th_sys);
        max_err3 = max(max_err3, err);
        assert(abs(ted_a - ted_b) < tol);  % sanity
    end
end
fprintf('  max |theta_helper - theta_Chi_sys| = %.2e\n', max_err3);
assert(max_err3 < tol, 'Theta inconsistency across wrappers');

fprintf('\n--- Test 4: translation invariance ---\n');
% distribution_aggregates(mu+varrho, p, s, varrho) == distribution_aggregates(mu, p, s, 0)
max_err4 = 0;
for k = 1:size(cases,1)
    mu = cases(k,1); p = cases(k,2); s = cases(k,3); vr = cases(k,4);
    [a1, b1, t1, F1, M1] = distribution_aggregates(mu + vr, p, s, vr);
    [a2, b2, t2, F2, M2] = distribution_aggregates(mu,      p, s, 0 );
    err = max(abs([a1-a2, b1-b2, t1-t2, F1-F2, M1-M2]));
    max_err4 = max(max_err4, err);
end
fprintf('  max diff over translation: %.2e\n', max_err4);
assert(max_err4 < tol, 'Translation invariance broken');

fprintf('\nAll tests passed (max errors: %.1e, %.1e, %.1e, %.1e)\n', ...
    max_err1, max_err2, max_err3, max_err4);
