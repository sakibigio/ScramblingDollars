%% test_sigma_grid.m
% Standalone check of the sigma grid / invariant distribution used by main_LFX,
% without running the whole model.  Sweeps (m_r1, m_r2, N) and reports whether
% the invariant distribution is well-behaved and whether the right tail is
% truncated relative to the FILTERED sigma range.

% Filtered sigma range under the tv-iota calibration (config C)
SIG_LO = 0.184; SIG_HI = 3.130;

pf  = fullfile('..','Filtering','data','MS_sigma_us_params.csv');
T   = readtable(pf);
gv  = @(nm) T.value(strcmp(T.param, nm));
mu1 = log(gv('sigma_ss_nor'));  rho1 = gv('rho_nor');  S1 = gv('Sigma_nor');
mu2 = log(gv('sigma_ss_scr'));  rho2 = gv('rho_scr');  S2 = gv('Sigma_scr');
p12 = gv('trans_nor'); p21 = gv('trans_scr');
P   = [1-p12 p12; p21 1-p21];
fprintf('Markov (new): mu1=%.4f rho1=%.4f S1=%.4f | mu2=%.4f rho2=%.4f S2=%.4f\n', ...
    mu1,rho1,S1,mu2,rho2,S2);
fprintf('uncond std: normal=%.4f  scrambling=%.4f\n', S1/sqrt(1-rho1^2), S2/sqrt(1-rho2^2));
fprintf('target: grid must span filtered sigma [%.3f, %.3f]\n\n', SIG_LO, SIG_HI);

specs = [ 8.0 2.5 102;   % current
          8.0 4.0 102;
          8.0 4.0 152;
          8.0 4.5 152;
          9.0 4.0 152;
          8.0 5.0 202 ];

fprintf('%5s %5s %5s | %16s %16s | %10s %10s | %8s %s\n', ...
    'm1','m2','N','normal sigma rng','scrambl sigma rng','endpt r1','endpt r2','valid','verdict');
fprintf('%s\n', repmat('-',1,118));
for k = 1:size(specs,1)
    m1 = specs(k,1); m2 = specs(k,2); N = specs(k,3);
    [Z1,Pi1] = tauchen(N/2, mu1, rho1, S1, m1);
    [Z2,Pi2] = tauchen(N/2, mu2, rho2, S2, m2);
    i12 = zeros(N/2,1); i21 = zeros(N/2,1);
    for ii = 1:N/2
        [~,i12(ii)] = min(abs(Z2-Z1(ii)));
        [~,i21(ii)] = min(abs(Z1-Z2(ii)));
    end
    P_full = zeros(N);
    for ii = 1:N/2
        P_full(ii,:)      = [P(1,1)*Pi1(ii,:)       P(1,2)*Pi2(i12(ii),:)];
        P_full(N/2+ii,:)  = [P(2,1)*Pi1(i21(ii),:)  P(2,2)*Pi2(ii,:)];
    end
    % ---- ROBUST invariant distribution: eigenvector for eigenvalue nearest 1
    [V,D] = eig(P_full');
    [~,j] = min(abs(diag(D)-1));
    v = real(V(:,j));
    if sum(v) < 0, v = -v; end
    v = max(v,0);  v = v/sum(v);
    v1 = v(1:N/2);      v1n = v1/sum(v1);
    v2 = v(N/2+1:end);  v2n = v2/sum(v2);
    valid = all(isfinite(v)) && abs(sum(v)-1) < 1e-8 && all(v >= -1e-12);
    lo1=exp(min(Z1)); hi1=exp(max(Z1)); lo2=exp(min(Z2)); hi2=exp(max(Z2));
    spans = (max(hi1,hi2) >= SIG_HI) && (min(lo1,lo2) <= SIG_LO);
    endpt1 = max(v1n(1), v1n(end));
    endpt2 = max(v2n(1), v2n(end));
    trunc  = (endpt1 > 1e-4) || (endpt2 > 1e-4);
    if ~valid,        verdict = 'INVALID invariant';
    elseif ~spans,    verdict = 'does not span filtered range';
    elseif trunc,     verdict = 'TRUNCATED (endpoint mass)';
    else,             verdict = 'OK';
    end
    fprintf('%5.1f %5.1f %5d | [%6.3f,%7.3f] [%6.3f,%7.3f] | %10.2e %10.2e | %8d %s\n', ...
        m1,m2,N, lo1,hi1, lo2,hi2, endpt1, endpt2, valid, verdict);
end
fprintf('\n(endpoint mass = max invariant mass on either extreme grid point, per regime)\n');
