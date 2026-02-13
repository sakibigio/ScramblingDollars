%% Test Unified Chi Functions
% Verifies that the unified functions with matching_type flag work correctly
clear; clc;

%% Parameters
theta_test = [0.3, 0.5, 0.8, 0.95, 1.0, 1.05, 1.5, 2.0];
lambda = 1.2;
eta = 0.5;
iota = 0.0012;  % 12 bps

%% Test Leontief (matching_type = 0)
fprintf('=== LEONTIEF (matching_type = 0) ===\n\n');
fprintf('theta\t\tChi_m\t\tChi_p\t\tPsi_m\t\tPsi_p\n');
fprintf('-----\t\t-----\t\t-----\t\t-----\t\t-----\n');

for th = theta_test
    cm = Chi_m(th, iota, lambda, eta, 0);
    cp = Chi_p(th, iota, lambda, eta, 0);
    pm = Psi_m(th, lambda, 0);
    pp = Psi_p(th, lambda, 0);
    fprintf('%.2f\t\t%.6f\t%.6f\t%.4f\t\t%.4f\n', th, cm, cp, pm, pp);
end

%% Test Cobb-Douglas (matching_type = 1)
fprintf('\n=== COBB-DOUGLAS (matching_type = 1) ===\n\n');
fprintf('theta\t\tChi_m\t\tChi_p\t\tPsi_m\t\tPsi_p\n');
fprintf('-----\t\t-----\t\t-----\t\t-----\t\t-----\n');

for th = theta_test
    cm = Chi_m(th, iota, lambda, eta, 1);
    cp = Chi_p(th, iota, lambda, eta, 1);
    pm = Psi_m(th, lambda, 1);
    pp = Psi_p(th, lambda, 1);
    fprintf('%.2f\t\t%.6f\t%.6f\t%.4f\t\t%.4f\n', th, cm, cp, pm, pp);
end

%% Test defaults (should be Leontief with eta=0.5)
fprintf('\n=== Testing Defaults ===\n');
th = 0.8;
cm_default = Chi_m(th, iota, lambda);  % Should use eta=0.5, Leontief
cm_explicit = Chi_m(th, iota, lambda, 0.5, 0);
fprintf('Default call: %.8f\n', cm_default);
fprintf('Explicit call: %.8f\n', cm_explicit);
fprintf('Match: %s\n', string(abs(cm_default - cm_explicit) < 1e-12));

%% Compare with old simplified formula (eta=0.5 hardcoded)
fprintf('\n=== Compare with Old Simplified Formula ===\n');
% Old formula from filter (eta=0.5 hardcoded, Leontief only)
chi_m_old = @(theta,iota,lambda) iota.*((theta+(1-theta).*exp(lambda)).^(1/2)-theta)./((1-theta).*exp(lambda));
chi_p_old = @(theta,iota,lambda) iota.*(theta.*(theta+(1-theta).*exp(lambda)).^(1/2)-theta)./((1-theta).*exp(lambda));

fprintf('theta\t\tChi_m (new)\tChi_m (old)\tMatch?\n');
for th = [0.3, 0.5, 0.8]
    cm_new = Chi_m(th, iota, lambda, 0.5, 0);
    cm_old = chi_m_old(th, iota, lambda);
    match = abs(cm_new - cm_old) < 1e-10;
    fprintf('%.2f\t\t%.8f\t%.8f\t%s\n', th, cm_new, cm_old, string(match));
end

fprintf('\nDone.\n');
