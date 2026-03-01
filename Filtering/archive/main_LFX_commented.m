% Archived commented-out code from main_LFX.m (lines 1110-1787)
% These were inactive plot experiments from earlier development.
% Archived on 2026-02-28 during code cleanup.
% ========================================================
%%%%
% cc=cc+1; plottype='_simexchguipRm';
% figure(cc)
% yyaxis left;plot(e_euus_t(500:end),'LineWidth',2);hold on;
% yyaxis right;plot(uip_Rm_t(500:end),'LineWidth',2);
% xlabel('Period','interpreter','latex');
% h=legend('Exchange Rate','$R^{eu}_{m}-R^{us}_{m}$ UIP Rm deviation','location','best');
% set(h,'interpreter','latex');
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 19, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Exchange Rate and UIP Rm deviation','interpreter','latex');
% if printit==1
%     orient landscape;
%     printsb(['fig' nameplot plottype]);
% end
% 
% cc=cc+1; plottype='_simuiplibor';
% figure(cc)
% plot(uip_dep_t(500:end),'LineWidth',2);xlabel('Period','interpreter','latex');
% ylabel('$R^{eu}_{libor}-R^{us}_{libor}$','interpreter','latex');
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 19, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('$R^{eu}_{libor}-R^{us}_{libor}$ (UIP libor deviation)','interpreter','latex');
% if printit==1
%     orient landscape;
%     printsb(['fig' nameplot plottype]);
% end

% %% Plot figures together
% % policy functions
% cc = 0;
% cc = cc+1; plottype='_policyfun';
% x_vec = sigma_us_vec(:);
% x_lab = '$\sigma_{us}$';
% % x_vec = (im_us_vec.^(freq)-1)*100;
% % x_lab = '$i^{us}_m(\%)$';
% figure(cc)
% subplot(4,4,1);
% splot1(x_vec,mu_us_vec);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\mu_{us}$','interpreter','latex');
% scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
% ylim(mean(mu_us_vec)+[-scale,scale]);
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Liquidity ratio (US)','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,2);
% splot1(x_vec,mu_eu_vec);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\mu_{eu}$','interpreter','latex');
% scale = max(0.001,max(abs(mu_eu_vec-mean(mu_eu_vec))));
% ylim(mean(mu_eu_vec)+[-scale,scale]);
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Liquidity ratio (EU)','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,3);
% splot1(x_vec,e_euus_vec);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('Exchange rate','interpreter','latex');
% scale = max(0.01,max(abs(e_euus_vec-mean(e_euus_vec))));
% ylim([mean(e_euus_vec)-scale,mean(e_euus_vec)+scale]);
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Exchange rate','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,4);
% splot1(x_vec,(Rm_us_vec.^(freq)-1)*100);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% ylim(([min([Rm_us_vec(:).^(freq);Rm_eu_vec(:).^(freq)]),max([Rm_us_vec(:).^(freq);Rm_eu_vec(:).^(freq)])]-1)*100);
% %scale = 10;
% %ylim(mean(([Rm_us_vec(:).^(freq);Rm_eu_vec(:).^(freq)]-1)*100)+[-scale,+scale]);
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_m$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,5);
% splot1(x_vec,(Rm_eu_vec.^(freq)-1)*100);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% ylim(([min([Rm_us_vec(:).^(freq);Rm_eu_vec(:).^(freq)]),max([Rm_us_vec(:).^(freq);Rm_eu_vec(:).^(freq)])]-1)*100);
% %scale = 0.001;
% %ylim(mean(([Rm_us_vec(:).^(freq);Rm_eu_vec(:).^(freq)]-1)*100)+[-scale,+scale]);
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_m$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,6);
% splot1(x_vec,(Rd_us_vec.^(freq)-1)*100);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% ylim(([min([Rd_us_vec(:).^(freq);Rd_eu_vec(:).^(freq)]),max([Rd_us_vec(:).^(freq);Rd_eu_vec(:).^(freq)])]-1)*100);
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_d$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,7);
% splot1(x_vec,(Rd_eu_vec.^(freq)-1)*100);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% ylim(([min([Rd_us_vec(:).^(freq);Rd_eu_vec(:).^(freq)]),max([Rd_us_vec(:).^(freq);Rd_eu_vec(:).^(freq)])]-1)*100);
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_d$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,8);
% splot1(x_vec,(RLibor_us_vec.^(freq)-1)*100);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% ylim(([min([RLibor_us_vec(:).^(freq);RLibor_eu_vec(:).^(freq)]),max([RLibor_us_vec(:).^(freq);RLibor_eu_vec(:).^(freq)])]-1)*100);
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_{libor}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,9);
% splot1(x_vec,(RLibor_eu_vec.^(freq)-1)*100);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% ylim(([min([RLibor_us_vec(:).^(freq);RLibor_eu_vec(:).^(freq)]),max([RLibor_us_vec(:).^(freq);RLibor_eu_vec(:).^(freq)])]-1)*100);
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{libor}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,10);
% splot1(x_vec,(RBond_us_vec.^(freq)-1)*100);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% ylim(([min([RBond_us_vec(:).^(freq);RBond_eu_vec(:).^(freq)]),max([RBond_us_vec(:).^(freq);RBond_eu_vec(:).^(freq)])]-1)*100);
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_{gov}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,11);
% splot1(x_vec,(RBond_eu_vec.^(freq)-1)*100);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% ylim(([min([RBond_us_vec(:).^(freq);RBond_eu_vec(:).^(freq)]),max([RBond_us_vec(:).^(freq);RBond_eu_vec(:).^(freq)])]-1)*100);
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{gov}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,12);
% splot1(x_vec,(Rb_us_vec.^(freq)-1)*100);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% ylim(([min([Rb_us_vec(:).^(freq);Rb_eu_vec(:).^(freq)]),max([Rb_us_vec(:).^(freq);Rb_eu_vec(:).^(freq)])]-1)*100);
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_b$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,13);
% splot1(x_vec,(Rb_eu_vec.^(freq)-1)*100);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% ylim(([min([Rb_us_vec(:).^(freq);Rb_eu_vec(:).^(freq)]),max([Rb_us_vec(:).^(freq);Rb_eu_vec(:).^(freq)])]-1)*100);
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_b$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,14);
% splot1(x_vec,(Rm_eu_vec.^(freq)-Rm_us_vec.^(freq))*100);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% scale=max(0.001,max(abs((Rm_eu_vec.^(freq)-Rm_us_vec.^(freq))*100-mean((Rm_eu_vec.^(freq)-Rm_us_vec.^(freq))*100))));
% %ylim(([min([Rm_eu_vec-Rm_us_vec]),max([Rm_eu_vec-Rm_us_vec])])*100);
% ylim([mean((Rm_eu_vec.^(freq)-Rm_us_vec.^(freq))*100)-scale,mean((Rm_eu_vec.^(freq)-Rm_us_vec.^(freq))*100)+scale]);
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{m}-R^{us}_{m}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,15);
% splot1(x_vec,(Rd_eu_vec.^(freq)-Rd_us_vec.^(freq))*100);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% %ylim(([min([Rd_eu_vec-Rd_us_vec]),max([Rd_eu_vec-Rd_us_vec])])*100);
% scale=max(0.001,max(abs((Rd_eu_vec.^(freq)-Rd_us_vec.^(freq))*100-mean((Rd_eu_vec.^(freq)-Rd_us_vec.^(freq))*100))));
% ylim([mean((Rd_eu_vec.^(freq)-Rd_us_vec.^(freq))*100)-scale,mean((Rd_eu_vec.^(freq)-Rd_us_vec.^(freq))*100)+scale]);
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{d}-R^{us}_{d}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,16);
% splot1(x_vec,(RLibor_eu_vec.^(freq)-RLibor_us_vec.^(freq))*100);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% %ylim(([min([RLibor_eu_vec-RLibor_us_vec]),max([RLibor_eu_vec-RLibor_us_vec])])*100);
% scale=max(0.001,max(abs((RLibor_eu_vec.^(freq)-RLibor_us_vec.^(freq))*100-mean((RLibor_eu_vec.^(freq)-RLibor_us_vec.^(freq))*100))));
% ylim([mean((RLibor_eu_vec.^(freq)-RLibor_us_vec.^(freq))*100)-scale,mean((RLibor_eu_vec.^(freq)-RLibor_us_vec.^(freq))*100)+scale]);
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{libor}-R^{us}_{libor}$','interpreter','latex','Fontsize',15);
% 
% if printit==1
%     orient landscape;
%     set(gcf, 'PaperUnits', 'normalized');
%     set(gcf, 'PaperPosition', [-0.30 -0.30 1.60 1.60]);
%     printsb(['fig' nameplot plottype]);
% end
% 
% % Quantities
% cc = cc+1;
% figure(cc); plottype='_policyquant';
% subplot(2,2,1);
% splot1(x_vec,d_us_vec);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('US real deposits $d_{us}$','interpreter','latex','Fontsize',15);
% 
% subplot(2,2,2);
% splot1(x_vec,d_eu_vec);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('EU real deposits $d_{eu}$','interpreter','latex','Fontsize',15);
% 
% subplot(2,2,3);
% splot1(x_vec,b_vec);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Real loans $b$','interpreter','latex','Fontsize',15);
% 
% if printit==1
%     orient landscape;
%     set(gcf, 'PaperUnits', 'normalized');
%     set(gcf, 'PaperPosition', [-0.06 -0.06 1.12 1.12]);    
%     printsb(['fig' nameplot plottype]);
% end
% 
% % Simulations
% cc=cc+1; plottype='_simpaths';
% Periods = 1:200;
% figure(cc)
% subplot(4,4,1);
% plot(Periods,mu_us_vec(chain(Periods)),'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\mu_{us}$','interpreter','latex');
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Liquidity ratio (US)','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,2);
% plot(Periods,mu_eu_vec(chain(Periods)),'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\mu_{eu}$','interpreter','latex');
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Liquidity ratio (EU)','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,3);
% plot(Periods,e_euus_vec(chain(Periods)),'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('Exchange rate','interpreter','latex');
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Exchange rate','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,4);
% plot(Periods,(Rm_us_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_m$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,5);
% plot(Periods,(Rm_eu_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_m$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,6);
% plot(Periods,(Rd_us_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_d$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,7);
% plot(Periods,(Rd_eu_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_d$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,8);
% plot(Periods,(RLibor_us_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_{libor}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,9);
% plot(Periods,(RLibor_eu_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{libor}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,10);
% plot(Periods,(RBond_us_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_{gov}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,11);
% plot(Periods,(RBond_eu_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{gov}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,12);
% plot(Periods,(Rb_us_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_b$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,13);
% plot(Periods,(Rb_eu_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_b$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,14);
% plot(Periods,(Rm_eu_vec(chain(Periods)).^(freq)-Rm_us_vec(chain(Periods)).^(freq))*100,'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{m}-R^{us}_{m}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,15);
% plot(Periods,(Rd_eu_vec(chain(Periods)).^(freq)-Rd_us_vec(chain(Periods)).^(freq))*100,'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{d}-R^{us}_{d}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,16);
% plot(Periods,(RLibor_eu_vec(chain(Periods)).^(freq)-RLibor_us_vec(chain(Periods)).^(freq))*100,'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{libor}-R^{us}_{libor}$','interpreter','latex','Fontsize',15);
% 
% if printit==1
%     orient landscape;
%     set(gcf, 'PaperUnits', 'normalized');
%     set(gcf, 'PaperPosition', [-0.06 -0.06 1.12 1.12]);    
%     printsb(['fig' nameplot plottype]);
% end

% sim paths of quantities
% cc=cc+1;plottype='_simpathsquant0';
% figure(cc);
% subplot(2,2,1);
% plot(Periods,d_us_vec(chain(Periods)),'Linewidth',2);
% ylabel('$d_{us}$','interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('US real deposits $d_{us}$','interpreter','latex','Fontsize',15);
% 
% subplot(2,2,2);
% plot(Periods,d_eu_vec(chain(Periods)),'Linewidth',2);
% ylabel('$d_{eu}$','interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('EU real deposits $d_{eu}$','interpreter','latex','Fontsize',15);
% 
% subplot(2,2,3);
% plot(Periods,b_vec(chain(Periods)),'Linewidth',2);
% ylabel('$b$','interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Real loans $b$','interpreter','latex','Fontsize',15);
% 
% if printit==1
%     orient landscape;
%     set(gcf, 'PaperUnits', 'normalized');
%     set(gcf, 'PaperPosition', [-0.06 -0.06 1.12 1.12]);    
%     printsb(['fig' nameplot plottype]);
% end
% 
% % Simulated shock paths
% 
% 
% cc=cc+1; plottype='_shockpath';
% figure(cc)
% subplot(1,2,1);
% plot(Periods,sigma_us_vec(chain(Periods)),'Linewidth',2);
% ylabel('$\sigma_{us}$','interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Dollar Volatility $\sigma_{us}$','interpreter','latex','Fontsize',15);
% 
% subplot(1,2,2);
% plot(Periods,(im_us_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);
% ylabel('$i^{us}_m(\%)$','interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Dollar Libor Rate $i^{us}_m(\%)$','interpreter','latex','Fontsize',15);
% 
% if printit==1
%     orient landscape;
%     set(gcf, 'Units','Inches','PaperUnits','Inches');
%     set(gcf, 'PaperSize', [11,4.5]);
%     set(gcf, 'PaperPosition', [0 0.25 11 4]);
%     printsb(['fig' nameplot plottype]);
%     print(gcf,'-dpdf',[path_g,filesep,['fig' nameplot plottype]])
% end
% 
% % Simulation 2
% % Simulations
% cc=cc+1; plottype='_simpaths2';
% Periods = 1:200;
% figure(cc)
% subplot(4,4,1);
% yyaxis left;plot(Periods,mu_us_vec(chain(Periods)),'Linewidth',2);
% ylabel('$\mu_{us}$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Liquidity ratio (US)','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,2);
% yyaxis left;plot(Periods,mu_eu_vec(chain(Periods)),'Linewidth',2);
% ylabel('$\mu_{eu}$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Liquidity ratio (EU)','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,3);
% yyaxis left;plot(Periods,e_euus_vec(chain(Periods)),'Linewidth',2);
% ylabel('Exchange rate','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Exchange rate','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,4);
% yyaxis left;plot(Periods,(Rm_us_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);
% ylabel('$\%$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_m$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,5);
% yyaxis left;plot(Periods,(Rm_eu_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);
% ylabel('$\%$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_m$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,6);
% yyaxis left;plot(Periods,(Rd_us_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);
% ylabel('$\%$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_d$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,7);
% yyaxis left;plot(Periods,(Rd_eu_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);
% ylabel('$\%$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_d$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,8);
% yyaxis left;plot(Periods,(RLibor_us_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);
% ylabel('$\%$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_{libor}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,9);
% yyaxis left;plot(Periods,(RLibor_eu_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);
% ylabel('$\%$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{libor}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,10);
% yyaxis left;plot(Periods,(RBond_us_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);
% ylabel('$\%$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_{gov}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,11);
% yyaxis left;plot(Periods,(RBond_eu_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);
% ylabel('$\%$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{gov}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,12);
% yyaxis left;plot(Periods,(Rb_us_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);
% ylabel('$\%$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_b$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,13);
% yyaxis left;plot(Periods,(Rb_eu_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);
% ylabel('$\%$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_b$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,14);
% yyaxis left;plot(Periods,(Rm_eu_vec(chain(Periods)).^(freq)-Rm_us_vec(chain(Periods)).^(freq))*100,'Linewidth',2);
% ylabel('$\%$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{m}-R^{us}_{m}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,15);
% yyaxis left;plot(Periods,(Rd_eu_vec(chain(Periods)).^(freq)-Rd_us_vec(chain(Periods)).^(freq))*100,'Linewidth',2);
% ylabel('$\%$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{d}-R^{us}_{d}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,16);
% yyaxis left;plot(Periods,(RLibor_eu_vec(chain(Periods)).^(freq)-RLibor_us_vec(chain(Periods)).^(freq))*100,'Linewidth',2);
% ylabel('$\%$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{libor}-R^{us}_{libor}$','interpreter','latex','Fontsize',15);
% 
% if printit==1
%     orient landscape;
%     set(gcf, 'PaperUnits', 'normalized');
%     set(gcf, 'PaperPosition', [-0.06 -0.06 1.12 1.12]);    
%     printsb(['fig' nameplot plottype]);
% end
% 
% % sim paths of quantities
% cc=cc+1;plottype='_simpathsquant';
% figure(cc);
% subplot(2,2,1);
% yyaxis left;plot(Periods,d_us_vec(chain(Periods)),'Linewidth',2);
% ylabel('$d_{us}$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('US real deposits $d_{us}$','interpreter','latex','Fontsize',15);
% 
% subplot(2,2,2);
% yyaxis left;plot(Periods,d_eu_vec(chain(Periods)),'Linewidth',2);
% ylabel('$d_{eu}$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('EU real deposits $d_{eu}$','interpreter','latex','Fontsize',15);
% 
% subplot(2,2,3);
% yyaxis left;plot(Periods,b_vec(chain(Periods)),'Linewidth',2);
% ylabel('$b$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Real loans $b$','interpreter','latex','Fontsize',15);
% 
% if printit==1
%     orient landscape;
%     set(gcf, 'PaperUnits', 'normalized');
%     set(gcf, 'PaperPosition', [-0.06 -0.06 1.12 1.12]);    
%     printsb(['fig' nameplot plottype]);
% end

