% simulate_model.m
% Extracted from main_LFX.m (lines 867–1109) on 2026-02-28
% Called by: main_LFX.m
%
% Contents:
%   - Setup plot vectors (X_vec, Y_vec from sigma grids)
%   - Monte Carlo simulation from Markov chain
%   - Correlations (exchange rate vs UIP deviations)
%   - Exchange rate regression (Engel-West style)
%   - Ergodic distribution and expected spreads

close all;
cc=0;

X_vec=sigma_eu_vec;
Y_vec=sigma_us_vec;

% NX=length(s_vec);
% NY=length(s_vec);
% X_mat=reshape(X_vec,NY,NX);
% Y_mat=reshape(Y_vec,NY,NX);
Z_vec=p_us_vec;
% Z_mat=reshape(Z_vec,NY,NX);

% cc=cc+1; plottype='_price_us';
% Z_vec=p_us_vec;
% figure(cc)
% splot1(Y_vec,Z_vec);
% %surf(X_mat,Y_mat,Z_mat);
% %xlabel('$\sigma_{eu}$','interpreter','latex');
% xlabel('$\sigma_{us}$','interpreter','latex');
% ylabel('$p_us$','interpreter','latex');
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 19, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Price Level U.S.','interpreter','latex');
% if printit==1
%     orient landscape;
%     printsb(['fig' nameplot plottype]);
% end
% 
% cc=cc+1; plottype='_exchangerate';
% Z_vec=e_euus_vec;
% %Z_mat=reshape(Z_vec,NY,NX);
% figure(cc)
% %surf(X_mat,Y_mat,Z_mat);
% splot1(Y_vec,Z_vec);
% %xlabel('$\sigma_{eu}$','interpreter','latex');
% xlabel('$\sigma_{us}$','interpreter','latex');
% ylabel('Exchange rate','interpreter','latex');
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 19, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Exchange rate','interpreter','latex');
% if printit==1
%     orient landscape;
%     printsb(['fig' nameplot plottype]);
% end
% 
% %cc=cc+1; plottype='_uip_dep';
% cc=cc+1; plottype='_uip';
% figure(cc)
% subplot(2,2,1);
% Z_vec=Rd_eu_vec-Rd_us_vec;
% %Z_mat=reshape(Z_vec,NY,NX);
% %surf(X_mat,Y_mat,Z_mat);
% splot1(Y_vec,Z_vec);
% %xlabel('$\sigma_{eu}$','interpreter','latex');
% xlabel('$\sigma_{us}$','interpreter','latex');
% ylabel('$R^{eu}_{d}-R^{us}_{d}$','interpreter','latex');
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 19, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('$R^{eu}_{d}-R^{us}_{d}$ (UIP deposit deviation)','interpreter','latex');
% % if printit==1
% %     orient landscape;
% %     printsb(['fig' nameplot plottype]);
% % end
% 
% % cc=cc+1; plottype='_uip_loan';
% subplot(2,2,2);
% Z_vec=Rb_eu_vec-Rb_us_vec;
% %Z_mat=reshape(Z_vec,NY,NX);
% % figure(cc)
% %surf(X_mat,Y_mat,Z_mat);
% splot1(Y_vec,Z_vec);
% %xlabel('$\sigma_{eu}$','interpreter','latex');
% xlabel('$\sigma_{us}$','interpreter','latex');
% ylabel('$R^{eu}_{b}-R^{us}_{b}$','interpreter','latex');
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 19, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('$R^{eu}_{b}-R^{us}_{b}$ (UIP loans deviation)','interpreter','latex');
% % if printit==1
% %     orient landscape;
% %     printsb(['fig' nameplot plottype]);
% % end
% 
% % cc=cc+1; plottype='_uip_bond';
% subplot(2,2,3);
% Z_vec=RBond_eu_vec-RBond_us_vec;
% %Z_mat=reshape(Z_vec,NY,NX);
% % figure(cc)
% %surf(X_mat,Y_mat,Z_mat);
% splot1(Y_vec,Z_vec);
% %xlabel('$\sigma_{eu}$','interpreter','latex');
% xlabel('$\sigma_{us}$','interpreter','latex');
% ylabel('$R^{eu}_{gov}-R^{us}_{gov}$','interpreter','latex');
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 19, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('$R^{eu}_{gov}-R^{us}_{gov}$ (UIP gov bond deviation)','interpreter','latex');
% % if printit==1
% %     orient landscape;
% %     printsb(['fig' nameplot plottype]);
% % end
% 
% % cc=cc+1; plottype='_uip_libor';
% subplot(2,2,4);
% Z_vec=RLibor_eu_vec-RLibor_us_vec;
% %Z_mat=reshape(Z_vec,NY,NX);
% % figure(cc)
% %surf(X_mat,Y_mat,Z_mat);
% splot1(Y_vec,Z_vec);
% %xlabel('$\sigma_{eu}$','interpreter','latex');
% xlabel('$\sigma_{us}$','interpreter','latex');
% ylabel('$R^{eu}_{libor}-R^{us}_{libor}$','interpreter','latex');
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 19, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('$R^{eu}_{libor}-R^{us}_{libor}$ (UIP libor deviation)','interpreter','latex');
% if printit==1
%     orient landscape;
%     printsb(['fig' nameplot plottype]);
% end
% 
% cc=cc+1; plottype='_uip_Rm';
% figure(cc)
% Z_vec=Rm_eu_vec-Rm_us_vec;
% splot1(Y_vec,Z_vec);
% xlabel('$\sigma_{us}$','interpreter','latex');
% ylabel('$R^{eu}_{m}-R^{us}_{m}$','interpreter','latex');
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 19, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('$R^{eu}_{m}-R^{us}_{m}$ (UIP Rm deviation)','interpreter','latex');
% if printit==1
%     orient landscape;
%     printsb(['fig' nameplot plottype]);
% end

% Simulate Markov Chain
rng('default');
start_value = round(size(Q_mat,1)/2);
chain_length = 5000;
chain = zeros(1,chain_length);
draw = zeros(1,chain_length);
chain(1) = start_value;

for t=2:chain_length
    this_step_dist = Q_mat(chain(t-1),:);
    cum_dist = cumsum(this_step_dist);
    draw(t) = rand();
    chain(t) = find(cum_dist>draw(t),1);
end


e_euus_t = e_euus_vec(chain(:));
bp_us_t  = Rb_us_vec(chain(:))-Rm_us_vec(chain(:));
bp_eu_t  = Rb_eu_vec(chain(:))-Rm_eu_vec(chain(:));
uip_Rm_t = Rm_eu_vec(chain(:))-Rm_us_vec(chain(:));
uip_dep_t = Rd_eu_vec(chain(:))-Rd_us_vec(chain(:));
uip_loan_t = Rb_eu_vec(chain(:))-Rb_us_vec(chain(:));
uip_bond_t = RBond_eu_vec(chain(:))-RBond_us_vec(chain(:));
uip_libor_t= RLibor_eu_vec(chain(:))-RLibor_us_vec(chain(:));
d_us_t = d_us_vec(chain(:));
d_eu_t = d_eu_vec(chain(:));
mu_us_t = mu_us_vec(chain(:));
mu_eu_t = mu_eu_vec(chain(:));
pidiff_t = pi_eu_vec(chain(:))-pi_us_vec(chain(:));
im_us_t = im_us_vec(chain(:));
im_eu_t = im_eu_vec(chain(:));
imdiff_t = im_eu_t-im_us_t;
pi_us_t  =pi_us_vec((chain(:)));

disp(['Corr. exchange rate vs UIP Rm: ' num2str(corr(e_euus_t',uip_Rm_t'))]);
disp(['Corr. exchange rate vs UIP Libor: ' num2str(corr(e_euus_t',uip_libor_t'))]);
disp(['Corr. exchange rate vs UIP Deposit: ' num2str(corr(e_euus_t',uip_dep_t'))]);
disp(['Corr. exchange rate vs UIP Loan: ' num2str(corr(e_euus_t',uip_loan_t'))]);
disp(['Corr. exchange rate vs UIP Bond: ' num2str(corr(e_euus_t',uip_bond_t'))]);



% Regression
%X_t = [diff(log(d_us_t(:))),diff(log(mu_us_t(:))),diff(uip_libor_t(:)),(pidiff_t(2:end)'),log(d_us_t(1:end-1)'),log(mu_us_t(1:end-1)),ones(999,1)];
%X_t = [diff(log(d_us_t(:))),diff(log(mu_us_t(:))),diff(uip_libor_t(:)),log(d_us_t(1:end-1)'),log(mu_us_t(1:end-1)),ones(999,1)];
% X_t = [diff(log(mu_us_t(:))),diff(uip_libor_t(:)),(pidiff_t(2:end)'),log(mu_us_t(1:end-1)'),ones(chain_length-1,1)];
%X_t = [diff(log(mu_us_t(:))),diff(uip_libor_t(:)),log(mu_us_t(1:end-1)'),ones(chain_length-1,1)];

%X_t = [diff(log(d_us_t(:))),diff(log(mu_us_t(:))),diff(imdiff_t(:)),(pidiff_t(2:end)'),log(d_us_t(1:end-1)'),log(mu_us_t(1:end-1)'),ones(chain_length-1,1)];
%X_t = [diff(log(d_us_t(:))),diff(log(mu_us_t(:))),diff(imdiff_t(:)),log(d_us_t(1:end-1)'),log(mu_us_t(1:end-1)'),ones(chain_length-1,1)];
%X_t = [diff(log(mu_us_t(:))),diff(imdiff_t(:)),(pidiff_t(2:end)'),log(mu_us_t(1:end-1)'),ones(chain_length-1,1)];
%X_t = [diff(log(mu_us_t(:))),diff(imdiff_t(:)),log(mu_us_t(1:end-1)'),ones(chain_length-1,1)];
%X_t = [diff(log(mu_us_t(:))),(pidiff_t(2:end)'),log(mu_us_t(1:end-1)'),ones(chain_length-1,1)];

%X_t = [diff(log(d_us_t(:))),diff(log(mu_us_t(:))),(pidiff_t(2:end)'),log(d_us_t(1:end-1)'),log(mu_us_t(1:end-1)'),ones(chain_length-1,1)];
%X_t = [diff(log(d_us_t(:))),diff(log(mu_us_t(:))),log(d_us_t(1:end-1)'),log(mu_us_t(1:end-1)'),ones(chain_length-1,1)];
X_t = [diff(log(mu_us_t(:))),(pidiff_t(2:end)'),log(mu_us_t(1:end-1)'),ones(chain_length-1,1)];
% X_t = [diff(log(mu_us_t(:))),log(mu_us_t(1:end-1)'),ones(chain_length-1,1)];

Y_t = diff(log(e_euus_t)');
sample = chain_length-191+1;
coef = (X_t(sample:end,:)'*X_t(sample:end,:))\(X_t(sample:end,:)'*Y_t(sample:end,:));
[coefb,cib,~,~,stats] = regress(Y_t(sample:end,:),X_t(sample:end,:));
talphaup = 1-0.05/2;
tnu = length(Y_t(sample:end,:))-length(coefb);
tvalue = tinv(talphaup,tnu);
sdb = (coefb-cib(:,1))/tvalue;
disp('         Coef   Std.Err.');
for i=1:length(coef)-1
    %disp(['beta',num2str(i),': ',num2str(coef(i)),',']);
    fprintf('beta%d: %0.4f  %0.4f\n',i,coefb(i),sdb(i));
end
%disp(['const: ',num2str(coef(end))]);
fprintf('const: %0.4f  %0.4f\n',coefb(end),sdb(end));
fprintf('R-sq: %0.4f\n',stats(1));

% moment table
ln_e_sample = log(e_euus_t(sample:end));
mean(ln_e_sample)
std(ln_e_sample)
corrcoef(ln_e_sample(2:end),ln_e_sample(1:end-1))
X=[ones(length(ln_e_sample)-1,1) ln_e_sample(1:end-1)'];
[B,BINT,R,RINT,STATS] = regress(ln_e_sample(2:end)' ,X);

disp('AR 1 coef of ln_e')
B(2)

disp('std residual')
sqrt( std(R)^2/(1-B(2)^2))

mean(log(d_us_t(sample:end)))
std(log(d_us_t(sample:end)))

mean(log(mu_us_t(sample:end)))
std(log(mu_us_t(sample:end)))
corrcoef(log(mu_us_t(sample:end-1))',log(mu_us_t(sample+1:end))')

mean(log(d_eu_t(sample:end)))
std(log(d_eu_t(sample:end)))

mean(log(mu_eu_t(sample:end)))
std(log(mu_eu_t(sample:end)))

mean(pidiff_t(sample:end))*100*12
std(pidiff_t(sample:end))*100*12
corrcoef(pidiff_t(sample:end-1)',pidiff_t(sample+1:end)')

b = zeros(size(Q_mat,1),1);
b(1) = 1;
row = [zeros(1,0),1,zeros(1,size(Q_mat,1)-1)];
A = eye(size(Q_mat))-Q_mat';
A(1,:) = row;
f_ss = A\b;
f_sum = sum(f_ss);
f_ss = f_ss./f_sum;
disp(['E[Rb_us-Rm_us] = ',num2str((Rb_us_vec.^(freq)-Rm_us_vec.^(freq))*f_ss(:)*1e4),' bps']);
disp(['E[Rm_eu-Rm_us] = ',num2str((Rm_eu_vec.^(freq)-Rm_us_vec.^(freq))*f_ss(:)*1e4),' bps']);
