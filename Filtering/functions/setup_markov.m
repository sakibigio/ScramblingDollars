% setup_markov.m
% Extracted from main_LFX.m (lines 337–476) on 2026-02-28
% Called by: main_LFX.m
%
% Contents:
%   - Markov regime parameters (normal + volatile)
%   - Tauchen discretization of sigma process
%   - Transition matrix construction (regime-switching)
%   - Invariant distribution computation
%   - State vector initialization
%
% NOTE: Update regime parameters below when re-estimating with CD matching.
%       Current values are from Leontief-filtered Markov estimation.

%% Global Solution
N_sigma_us = 102;

% Rates
Z_im_us = (imss_us);
Zprob_im_us = 1    ;
Z_im_eu = (imss_eu);
Zprob_im_eu = 1    ;

% Index by States
index1=1:N_sigma_us/2;
index2=N_sigma_us/2+1:N_sigma_us;

% % % Normal Regime:
mu_sigma_us_r1 = -1.27 ;
rho_sigma_us_r1 = 0.88 ;
Sigma_sigma_us_r1 = 0.2;
m_sigma_us_r1 = 8.0   ;
% 
% % Volatile Regime
mu_sigma_us_r2 = -0.799;
rho_sigma_us_r2 = 0.61 ;
Sigma_sigma_us_r2 = 0.9;
m_sigma_us_r2 = 2.5    ;

% % Transition Probabilities
P=[0.984 0.016; 0.061 0.939];

% % % Normal Regime:
% mu_sigma_us_r1 = -1.0;
% rho_sigma_us_r1 = 0.975;
% Sigma_sigma_us_r1 = 0.2;
% m_sigma_us_r1 = 3.5;
% 
% % Volatile Regime
% mu_sigma_us_r2 = -1.0;
% rho_sigma_us_r2 = 0.0;
% Sigma_sigma_us_r2 = 0.2;
% m_sigma_us_r2 = 14.5;
% 
% % Transition Probabilities
% P=[0.999 0.001; 0.001 0.999];

% Construction of process
Tauchen_method=1; Rouwenhorst_method=0;
if Tauchen_method==1
    [Z_sigma_us_r1,Zprob_sigma_us_r1] = tauchen(N_sigma_us/2,mu_sigma_us_r1,rho_sigma_us_r1,Sigma_sigma_us_r1,m_sigma_us_r1);
    [Z_sigma_us_r2,Zprob_sigma_us_r2] = tauchen(N_sigma_us/2,mu_sigma_us_r2,rho_sigma_us_r2,Sigma_sigma_us_r2,m_sigma_us_r2);

elseif Rouwenhorst_method==1
    [Z_sigma_us_r1,Zprob_sigma_us_r1] =discretizeAR1_Rouwenhorst(mu_sigma_us_r1*(1-rho_sigma_us_r1),rho_sigma_us_r1,Sigma_sigma_us_r1,N_sigma_us/2);
    [Z_sigma_us_r2,Zprob_sigma_us_r2] =discretizeAR1_Rouwenhorst(mu_sigma_us_r2,rho_sigma_us_r2,Sigma_sigma_us_r2,N_sigma_us/2);
end
[eig_vecs1,eigs1]=eig(Zprob_sigma_us_r1');
invp1=eig_vecs1(:,1)/(sum(eig_vecs1(:,1)));
[eig_vecs2,eigs2]=eig(Zprob_sigma_us_r2');
invp2=eig_vecs2(:,1)/(sum(eig_vecs2(:,1)));
figure('Name','Invariant Distributions')
area(exp(Z_sigma_us_r1),invp1,'FaceAlpha',0.5,'FaceColor','b'); hold on;
area(exp(Z_sigma_us_r2),invp2,'FaceAlpha',0.5,'FaceColor','r'); hold on;

% Test location
figure('Name','Locations')
scatter(Z_sigma_us_r1,Z_sigma_us_r1*0+1,'Filled','Color',[0.1 0.5 0.8]); hold on;
scatter(Z_sigma_us_r2,Z_sigma_us_r2*0+1.1,'Filled','Color',[0.8 0.2 0.8]); grid on;
legend('low','high');

% Index Key - Find nearest location for construction of process:
for ii=1:N_sigma_us/2
    [~,index_r1_to_r2(ii)]=min(abs(Z_sigma_us_r2-Z_sigma_us_r1(ii)));
    [~,index_r2_to_r1(ii)]=min(abs(Z_sigma_us_r1-Z_sigma_us_r2(ii)));
end

% Single Chain
Z_sigma_us=[Z_sigma_us_r1; Z_sigma_us_r2];
Zprob_sigma_us=zeros(N_sigma_us);
for ii=1:N_sigma_us/2
    Zprob_sigma_us(ii,:)=[P(1,1)*Zprob_sigma_us_r1(ii,:) P(1,2)*Zprob_sigma_us_r2(index_r1_to_r2(ii),:)];
    Zprob_sigma_us(N_sigma_us/2+ii,:)=[P(2,1)*Zprob_sigma_us_r1(index_r2_to_r1(ii),:) P(2,2)*Zprob_sigma_us_r2(ii,:)];
end
[eig_vecs,eigs]=eig(Zprob_sigma_us');
invp=eig_vecs(:,1)/(sum(eig_vecs(:,1)));
invpp1=invp(1:N_sigma_us/2); invp1=invpp1/sum(invpp1);
invpp2=invp(N_sigma_us/2+1:end); invp2=invpp2/sum(invpp2);
figure('Name','Conditional Distributions')
area(exp(Z_sigma_us_r1),invp1,'FaceAlpha',0.5,'FaceColor','b'); hold on;
area(exp(Z_sigma_us_r2),invp2,'FaceAlpha',0.5,'FaceColor','r'); hold on;

% Other Parameters
N_Theta_d_us = 1;
mu_Theta_d_us = 0;
%rho_Theta_d_us = 0.98;
rho_Theta_d_us = 0.9828;
%sigma_Theta_d_us = 0.05;
sigma_Theta_d_us = 0.073;
m_Theta_d_us = 3;
[Z_Theta_d_us,Zprob_Theta_d_us] = tauchen(N_Theta_d_us,mu_Theta_d_us,rho_Theta_d_us,sigma_Theta_d_us,m_Theta_d_us);
Z_Theta_d_us = 0;
Zprob_Theta_d_us = 1;

% N_sigma_eu = 11;
% mu_sigma_eu = log(sigma_eu);
% rho_sigma_eu = 0.95;
% sigma_sigma_eu = 0.2;
% m_sigma_eu = 3;
% [Z_sigma_eu,Zprob_sigma_eu] = tauchen(N_sigma_eu,mu_sigma_eu,rho_sigma_eu,sigma_sigma_eu,m_sigma_eu);

sigmashock_us_vec = zeros(1,length(Z_sigma_us)*length(Z_im_us)*length(Z_im_eu)*length(Z_Theta_d_us));
imshock_us_vec = zeros(1,length(Z_sigma_us)*length(Z_im_us)*length(Z_im_eu)*length(Z_Theta_d_us));
imshock_eu_vec = zeros(1,length(Z_sigma_us)*length(Z_im_us)*length(Z_im_eu)*length(Z_Theta_d_us));
Thetadshock_us_vec = zeros(1,length(Z_sigma_us)*length(Z_im_us)*length(Z_im_eu)*length(Z_Theta_d_us));
Q_mat = zeros(length(Z_sigma_us)*length(Z_im_us)*length(Z_im_eu)*length(Z_Theta_d_us),length(Z_sigma_us)*length(Z_im_us)*length(Z_im_eu)*length(Z_Theta_d_us));
id = 0;
for i1=1:length(Z_im_us)
    for i2=1:length(Z_sigma_us)
        for i3=1:length(Z_im_eu)
            for i4=1:length(Z_Theta_d_us)
                id = id+1;
                sigmashock_us_vec(id) = exp(Z_sigma_us(i2));
                imshock_us_vec(id) = Z_im_us(i1);
                imshock_eu_vec(id) = Z_im_eu(i3);
                Thetadshock_us_vec(id) = exp(Z_Theta_d_us(i4));
                temp = Zprob_sigma_us(i2,:)'*Zprob_im_us(i1,:);
                temp = Zprob_im_eu(i3,:)'*temp(:)';
                temp = Zprob_Theta_d_us(i4,:)'*temp(:)';
                Q_mat(id,:) = temp(:)';
            end
        end
    end
end

Q_mat(Q_mat<1e-9) = 0;
Q_mat = sparse(Q_mat);
N_s=size(Q_mat,1);

if N_s==length(Z_sigma_us)
    Q_mat=Zprob_sigma_us;
end
% Initiate Endogenous Variables
LFX_empty_vecs_state; % initiates times series
