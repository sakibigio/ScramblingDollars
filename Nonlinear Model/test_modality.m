%% test_modality.m - is the "normal regime" invariant distribution bimodal?
% An AR(1) invariant should be unimodal.  But invp1 as computed in main_LFX is
% the regime-1 SLICE OF THE JOINT CHAIN, which mixes (a) mass that has been in
% regime 1 a long time with (b) mass that just switched in from regime 2 at the
% nearest-neighbour mapped locations.  If the regime means differ, that can
% create a second mode.  This script checks, for a given params file.

if ~exist('pfile','var'), pfile = '../Filtering/data/MS_sigma_us_params_C.csv'; end
if ~exist('m1','var'), m1 = 8; end
if ~exist('m2','var'), m2 = 4.5; end
if ~exist('N','var'),  N  = 152; end

T  = readtable(pfile);
gv = @(nm) T.value(strcmp(T.param,nm));
mu1=log(gv('sigma_ss_nor')); rho1=gv('rho_nor'); S1=gv('Sigma_nor');
mu2=log(gv('sigma_ss_scr')); rho2=gv('rho_scr'); S2=gv('Sigma_scr');
p12=gv('trans_nor'); p21=gv('trans_scr'); P=[1-p12 p12; p21 1-p21];
fprintf('params: %s\n', pfile);
fprintf('  normal: mu=%.4f rho=%.4f S=%.4f | scrambling: mu=%.4f rho=%.4f S=%.4f\n', ...
    mu1,rho1,S1,mu2,rho2,S2);

[Z1,Pi1] = tauchen(N/2, mu1, rho1, S1, m1);
[Z2,Pi2] = tauchen(N/2, mu2, rho2, S2, m2);

% --- standalone AR(1) invariant for regime 1 (the theoretical benchmark) ---
[V1,D1]=eig(Pi1'); [~,j1]=min(abs(diag(D1)-1));
a1=real(V1(:,j1)); if sum(a1)<0, a1=-a1; end; a1=max(a1,0); a1=a1/sum(a1);

% --- joint chain, then slice regime 1 (what main_LFX actually plots) ---
i12=zeros(N/2,1); i21=zeros(N/2,1);
for ii=1:N/2
    [~,i12(ii)]=min(abs(Z2-Z1(ii)));
    [~,i21(ii)]=min(abs(Z1-Z2(ii)));
end
Pf=zeros(N);
for ii=1:N/2
    Pf(ii,:)     =[P(1,1)*Pi1(ii,:)      P(1,2)*Pi2(i12(ii),:)];
    Pf(N/2+ii,:) =[P(2,1)*Pi1(i21(ii),:) P(2,2)*Pi2(ii,:)];
end
[V,D]=eig(Pf'); [~,jj]=min(abs(diag(D)-1));
v=real(V(:,jj)); if sum(v)<0, v=-v; end; v=max(v,0); v=v/sum(v);
b1=v(1:N/2); b1=b1/sum(b1);       % regime-1 slice of the JOINT invariant
b2=v(N/2+1:end); b2=b2/sum(b2);

    function n = nmodes(p, tol)
        % count interior local maxima with a relative prominence filter
        n = 0; pk = [];
        for k = 2:numel(p)-1
            if p(k) > p(k-1) && p(k) >= p(k+1), pk(end+1) = k; end %#ok<AGROW>
        end
        % keep peaks that are at least tol * max(p)
        pk = pk(p(pk) > tol*max(p));
        n = numel(pk);
    end

tol = 0.02;
fprintf('\nmodes (peaks > %.0f%% of max):\n', tol*100);
fprintf('  standalone AR(1) invariant, NORMAL     : %d\n', nmodes(a1,tol));
fprintf('  joint-chain slice,          NORMAL     : %d   <-- what main_LFX plots\n', nmodes(b1,tol));
fprintf('  joint-chain slice,          SCRAMBLING : %d\n', nmodes(b2,tol));

% locate the peaks in sigma units
pk = [];
for k=2:numel(b1)-1
    if b1(k)>b1(k-1) && b1(k)>=b1(k+1) && b1(k)>tol*max(b1), pk(end+1)=k; end %#ok<AGROW>
end
if numel(pk)>1
    fprintf('\n  NORMAL-slice peak locations (sigma): %s\n', mat2str(round(exp(Z1(pk))',4)));
    fprintf('  peak heights (rel. to max): %s\n', mat2str(round(b1(pk)'/max(b1),3)));
    fprintf('  scrambling regime mean sigma = %.4f  (nearest r1 grid pt = %.4f)\n', ...
        exp(mu2), exp(Z1(min(abs(Z1-mu2))==abs(Z1-mu2))));
end
fprintf('\n  E[sigma]: standalone-normal=%.4f  joint-normal=%.4f  joint-scrambling=%.4f\n', ...
    exp(Z1(:))'*a1(:), exp(Z1(:))'*b1(:), exp(Z2(:))'*b2(:));
