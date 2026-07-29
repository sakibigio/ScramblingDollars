%% test_eig_live.m
% Replicates functions/setup_markov.m's chain construction (LIVE Cobb-Douglas
% pipeline) and checks whether its eig-ordering assumption holds.
% setup_markov lines 108/110/133 take eig_vecs(:,1) as the invariant
% distribution, but eig() does not order eigenvalues -- column 1 is correct
% only by luck.
addpath('functions'); addpath('utils'); addpath('data');

f = fullfile('data','MS_params.csv');
if ~isfile(f), error('MS_params.csv not found'); end
fid=fopen(f,'r'); fgetl(fid); d=textscan(fid,'%s%f','Delimiter',','); fclose(fid);
M = containers.Map(d{1}, num2cell(d{2}));
mu1=M('mu_sigma_us_r1'); rho1=M('rho_sigma_us_r1'); S1=M('Sigma_sigma_us_r1');
mu2=M('mu_sigma_us_r2'); rho2=M('rho_sigma_us_r2'); S2=M('Sigma_sigma_us_r2');
P=[M('P11') M('P12'); M('P21') M('P22')];
N_sigma_us = 152; m_common = 5.0;
fprintf('LIVE MS params: r1 mu=%.4f rho=%.4f S=%.4f | r2 mu=%.4f rho=%.4f S=%.4f\n', ...
    mu1,rho1,S1,mu2,rho2,S2);
fprintf('P = [%.4f %.4f; %.4f %.4f]\n', P(1,1),P(1,2),P(2,1),P(2,2));

su1=S1/sqrt(max(1e-6,1-rho1^2)); su2=S2/sqrt(max(1e-6,1-rho2^2));
mu_c=(mu1+mu2)/2; su_c=max(su1,su2);
[Zc,~]=tauchen(N_sigma_us/2, mu_c, 0, su_c, m_common);
Ng=N_sigma_us/2; a1=(1-rho1)*mu1; a2=(1-rho2)*mu2; zs=Zc(2)-Zc(1);
Pi1=zeros(Ng); Pi2=zeros(Ng);
for j=1:Ng
  for k=1:Ng
    if k==1
      Pi1(j,k)=0.5*erfc(-(Zc(1)-a1-rho1*Zc(j)+zs/2)/(S1*sqrt(2)));
      Pi2(j,k)=0.5*erfc(-(Zc(1)-a2-rho2*Zc(j)+zs/2)/(S2*sqrt(2)));
    elseif k==Ng
      Pi1(j,k)=1-0.5*erfc(-(Zc(Ng)-a1-rho1*Zc(j)-zs/2)/(S1*sqrt(2)));
      Pi2(j,k)=1-0.5*erfc(-(Zc(Ng)-a2-rho2*Zc(j)-zs/2)/(S2*sqrt(2)));
    else
      Pi1(j,k)=0.5*erfc(-(Zc(k)-a1-rho1*Zc(j)+zs/2)/(S1*sqrt(2)))-0.5*erfc(-(Zc(k)-a1-rho1*Zc(j)-zs/2)/(S1*sqrt(2)));
      Pi2(j,k)=0.5*erfc(-(Zc(k)-a2-rho2*Zc(j)+zs/2)/(S2*sqrt(2)))-0.5*erfc(-(Zc(k)-a2-rho2*Zc(j)-zs/2)/(S2*sqrt(2)));
    end
  end
end
Pf=zeros(N_sigma_us);
for ii=1:Ng
  Pf(ii,:)   =[P(1,1)*Pi1(ii,:) P(1,2)*Pi2(ii,:)];
  Pf(Ng+ii,:)=[P(2,1)*Pi1(ii,:) P(2,2)*Pi2(ii,:)];
end

fprintf('\n===== eig-ordering check (LIVE pipeline) =====\n');
names={'regime 1','regime 2','joint chain'}; mats={Pi1,Pi2,Pf};
for k=1:3
  [V,D]=eig(mats{k}'); dd=diag(D);
  [~,js]=min(abs(dd-1));
  v1=real(V(:,1)); s1=sum(v1);
  status='CORRECT'; if js~=1, status='*** WRONG ***'; end
  fprintf('%-12s: col1 eigval=%+.6f | unit-eigval in col %3d | %s\n', ...
      names{k}, real(dd(1)), js, status);
  if js~=1
    v1n=v1/s1;
    vs=real(V(:,js)); if sum(vs)<0, vs=-vs; end; vs=max(vs,0); vs=vs/sum(vs);
    fprintf('              col1 has %d negative entries; max|col1-true| = %.3e\n', ...
        sum(v1n<-1e-10), max(abs(v1n(:)-vs(:))));
  end
end

[V,D]=eig(Pf'); [~,js]=min(abs(diag(D)-1));
v=real(V(:,js)); if sum(v)<0, v=-v; end; v=max(v,0); v=v/sum(v);
A=v(1:Ng); A=A/sum(A); B=v(Ng+1:end); B=B/sum(B);
fprintf('\nGrid: N=%d, sigma range [%.4f, %.4f] (m_common=%.1f)\n', ...
    N_sigma_us, exp(min(Zc)), exp(max(Zc)), m_common);
fprintf('TRUE invariant endpoint mass: r1 lo=%.2e hi=%.2e | r2 lo=%.2e hi=%.2e\n', ...
    A(1),A(end),B(1),B(end));
nm=@(p) sum(p(2:end-1)>p(1:end-2) & p(2:end-1)>=p(3:end) & p(2:end-1)>0.02*max(p));
fprintf('modes (>2%% of max): regime1=%d  regime2=%d\n', nm(A), nm(B));
hl1=log(0.5)/log(rho1); dur1=1/P(1,2);
fprintf('rho1 half-life=%.1f mo vs regime-1 duration=%.1f mo -> %s\n', hl1, dur1, ...
    string(hl1>dur1));
