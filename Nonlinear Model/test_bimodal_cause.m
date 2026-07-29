%% test_bimodal_cause.m
% Is the bimodal normal-regime slice a coding bug, or a consequence of the
% estimated persistence?  Sweep rho_normal holding everything else at C's
% values and find where the second mode appears.
%
% Hypothesis: the joint chain carries the LEVEL across regime switches
% (index_r2_to_r1).  Mass entering regime 1 near the scrambling mean relaxes
% toward the normal mean at within-regime half-life  ln(.5)/ln(rho).  If that
% half-life exceeds the expected regime duration 1/trans_nor, the mass never
% relaxes and shows up as a separate mode.

T  = readtable('../Filtering/data/MS_sigma_us_params_C.csv');
gv = @(nm) T.value(strcmp(T.param,nm));
mu1=log(gv('sigma_ss_nor')); S1=gv('Sigma_nor');
mu2=log(gv('sigma_ss_scr')); rho2=gv('rho_scr'); S2=gv('Sigma_scr');
p12=gv('trans_nor'); p21=gv('trans_scr');
rho1_hat = gv('rho_nor'); rho1_se = gv('rho_nor_se');
N=152; m1=8; m2=4.5;

fprintf('C estimate: rho_normal = %.4f  (s.e. %.4f)  -> 95%% CI [%.4f, %.4f]\n', ...
    rho1_hat, rho1_se, rho1_hat-1.96*rho1_se, rho1_hat+1.96*rho1_se);
fprintf('expected regime-1 duration = 1/trans_nor = %.0f months\n\n', 1/p12);

fprintf('%8s %12s %12s %8s %s\n','rho_nor','half-life','duration','modes','');
fprintf('%s\n', repmat('-',1,58));
for rho1 = [0.90 0.95 0.97 0.9729 0.98 0.99 0.995 rho1_hat 0.999]
    P=[1-p12 p12; p21 1-p21];
    [Z1,Pi1]=tauchen(N/2,mu1,rho1,S1,m1);
    [Z2,Pi2]=tauchen(N/2,mu2,rho2,S2,m2);
    i12=zeros(N/2,1); i21=zeros(N/2,1);
    for ii=1:N/2
        [~,i12(ii)]=min(abs(Z2-Z1(ii)));
        [~,i21(ii)]=min(abs(Z1-Z2(ii)));
    end
    Pf=zeros(N);
    for ii=1:N/2
        Pf(ii,:)    =[P(1,1)*Pi1(ii,:)      P(1,2)*Pi2(i12(ii),:)];
        Pf(N/2+ii,:)=[P(2,1)*Pi1(i21(ii),:) P(2,2)*Pi2(ii,:)];
    end
    [V,D]=eig(Pf'); [~,jj]=min(abs(diag(D)-1));
    v=real(V(:,jj)); if sum(v)<0, v=-v; end; v=max(v,0); v=v/sum(v);
    b1=v(1:N/2); b1=b1/sum(b1);
    nm=0;
    for k=2:numel(b1)-1
        if b1(k)>b1(k-1) && b1(k)>=b1(k+1) && b1(k)>0.02*max(b1), nm=nm+1; end
    end
    hl = log(0.5)/log(rho1);
    mark = '';
    if abs(rho1-rho1_hat)<1e-9, mark=' <- C estimate'; end
    if abs(rho1-0.9729)<1e-9,   mark=' <- B estimate'; end
    fprintf('%8.4f %12.1f %12.0f %8d%s\n', rho1, hl, 1/p12, nm, mark);
end
fprintf('\n(half-life and duration in months; modes = peaks >2%% of max in the normal slice)\n');
