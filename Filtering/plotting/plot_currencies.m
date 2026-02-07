%% RW_filter_plot_baseline_curr
% plot baseline results from filtering excercise
% for other currencies
%% Plotting Other Currencies
figure('Name','Diagnostics','NumberTitle','off') 
subplot(3,3,1); 
plot(dates(datesperiod),sigma_eu_TED_flag,'LineWidth',2);hold on;
datetick('x','yy');
grid on; axis tight;
label_x('Time (Year)');
formataxis(gca);
h=legend('Country','US');
set(h,'interpreter','latex','location','Northeast','Fontsize',10, ...
    'box','off');
title('EU','interpreter','latex','Fontsize',15);
for j=1:length(curlist)
    aux = eval(['sigma_' curlist{j} '_TED_flag(datesperiod)']);
    subplot(3,3,j+1);
    plot(dates(datesperiod),aux,'LineWidth',2); hold on;
    %subplot(3,3,j+1);plot(eval(['-inv_e_' curlist{j} '(end-200:end)']),'LineWidth',2);
    datetick('x','yy');
    grid on; axis tight;
    label_x('Time (Year)');
    formataxis(gca);
    title([conlist{j}],'interpreter','latex','Fontsize',15);
end

% Sigma per currencie
figure('Name','Currency Sigmas','NumberTitle','off') 
subplot(3,3,1); 
plot(dates(datesperiod),log(sigma_eu_TED_t(datesperiod)),'LineWidth',2);hold on;
plot(dates(datesperiod),log(sigma_us_TED_t(datesperiod)),'LineWidth',2);hold on;
datetick('x','yy');
grid on; axis tight;
label_x('Time (Year)');
formataxis(gca);
h=legend('Country','US');
set(h,'interpreter','latex','location','Northeast','Fontsize',10,'box','off');
title('EU','interpreter','latex','Fontsize',15);
for j=1:length(curlist)
    aux = eval(['sigma_' curlist{j} '_t(datesperiod)']);
    %temp = eval(['ln_' curlist{j} '_us_t((datesperiod(1)))'])-eval(['-oo_.UpdatedVariables.inv_e_' curlist{j} '(datesperiod(1))']);
    subplot(3,3,j+1);
    plot(dates(datesperiod),log(aux(datesperiod)),'LineWidth',2); hold on;
    plot(dates(datesperiod),log(sigma_us_TED_t(datesperiod)),'LineWidth',1,'LineStyle','-');
    %subplot(3,3,j+1);plot(eval(['-inv_e_' curlist{j} '(end-200:end)']),'LineWidth',2);
    datetick('x','yy');
    grid on; axis tight;
    label_x('Time (Year)');
    formataxis(gca);
    title([conlist{j}],'interpreter','latex','Fontsize',15);
end

figure('Name','Global Bond Premia','NumberTitle','off') 
subplot(3,3,1); 
plot(dates(datesperiod),BP_eu_t(datesperiod)*abs_scale,'LineWidth',2);hold on;
plot(dates(datesperiod),BP_s_eu_t(datesperiod)*abs_scale,'LineWidth',2);hold on;
datetick('x','yy');
grid on; axis tight;
label_x('Time (Year)');
formataxis(gca);
h=legend('Model','Data');
set(h,'interpreter','latex','location','Northeast','Fontsize',10);
title('EU','interpreter','latex','Fontsize',15);
for j=1:length(curlist)
    aux = eval(['BP_' curlist{j} '_t(datesperiod)']);
    aux_d = eval(['BP_s_' curlist{j} '_t(datesperiod)']);
    %temp = eval(['ln_' curlist{j} '_us_t((datesperiod(1)))'])-eval(['-oo_.UpdatedVariables.inv_e_' curlist{j} '(datesperiod(1))']);
    subplot(3,3,j+1);
    plot(dates(datesperiod),aux(datesperiod)*abs_scale,'LineWidth',2); hold on;
    plot(dates(datesperiod),aux_d(datesperiod)*abs_scale,'LineWidth',2);
%    plot(dates(datesperiod),aux_d(datesperiod),'LineWidth',2);
    %subplot(3,3,j+1);plot(eval(['-inv_e_' curlist{j} '(end-200:end)']),'LineWidth',2);
    datetick('x','yy');
    grid on; axis tight;
    label_x('Time (Year)');
    formataxis(gca);
    title([conlist{j}],'interpreter','latex','Fontsize',15);
end

% CIP Deviations
figure('Name','CIP Deviations','NumberTitle','off') 
subplot(3,3,1); 
plot(dates(datesperiod),CIP_t(datesperiod)*abs_scale,'LineWidth',2);hold on;
plot(dates(datesperiod),CIP_s_eu_t(datesperiod)*abs_scale,'LineWidth',2);hold on;
datetick('x','yy');
grid on; axis tight;
label_x('Time (Year)');
formataxis(gca);
h=legend('Model','Data');
set(h,'interpreter','latex','location','Northeast','Fontsize',10);
title('EU','interpreter','latex','Fontsize',15);
for j=1:length(curlist)
    aux = eval(['CIP_' curlist{j} '_t(datesperiod)']);
    aux_d = eval(['CIP_s_' curlist{j} '_t(datesperiod)']);
    %temp = eval(['ln_' curlist{j} '_us_t((datesperiod(1)))'])-eval(['-oo_.UpdatedVariables.inv_e_' curlist{j} '(datesperiod(1))']);
    subplot(3,3,j+1);
    plot(dates(datesperiod),aux(datesperiod)*abs_scale,'LineWidth',2); hold on;
    plot(dates(datesperiod),aux_d(datesperiod)*abs_scale,'LineWidth',2);
    %subplot(3,3,j+1);plot(eval(['-inv_e_' curlist{j} '(end-200:end)']),'LineWidth',2);
    datetick('x','yy');
    grid on; axis tight;
    label_x('Time (Year)');
    formataxis(gca);
    title([conlist{j}],'interpreter','latex','Fontsize',15);
end

% FX deviations
figure('Name','FX Rates','NumberTitle','off')
temp = mean(ln_eu_us_t(datesperiod))+mean(log(inv_e_t(datesperiod)));
subplot(3,3,1);
plot(dates(datesperiod),-log(inv_e_t(datesperiod))+temp,'LineWidth',2);hold on;
plot(dates(datesperiod),ln_eu_us_t(datesperiod),'LineWidth',2,'LineStyle','--');hold off;
datetick('x','yy');
grid on; axis tight;
label_x('Time (Year)');
formataxis(gca);
h=legend('Model','Data');
set(h,'interpreter','latex','location','Northeast','Fontsize',10);
title('EUR/USD','interpreter','latex','Fontsize',15);
for j=1:length(curlist)
    temp = mean(eval(['ln_' curlist{j} '_us_t((datesperiod))']))+mean(eval(['log(inv_e_' curlist{j} '_t)']));
    subplot(3,3,j+1);
    plot(dates(datesperiod),eval(['-log(inv_e_' curlist{j} '_t)'])+temp,'LineWidth',2);hold on;
    plot(dates(datesperiod),eval(['ln_' curlist{j} '_us_t((datesperiod))']),'LineWidth',2);hold off;
    datetick('x','yy');
    grid on; axis tight;
    label_x('Time (Year)');
    formataxis(gca);
    title([conlist{j} '/USD'],'interpreter','latex','Fontsize',15);
end

%% Multicountry Plots
figure('Name','FX Rates','NumberTitle','off')
temp1 = mean(ln_eu_us_t(datesperiod));
temp2 = mean(log(inv_e_t(datesperiod)));
subplot(3,3,1);
plot(dates(datesperiod),-log(inv_e_t(datesperiod))+temp2,'LineWidth',2);hold on;
plot(dates(datesperiod),ln_eu_us_t(datesperiod)-temp1,'LineWidth',2,'LineStyle','--');hold off;
datetick('x','yy');
grid on; axis tight;
label_x('Time (Year)');
formataxis(gca);
h=legend('Model','Data');
set(h,'interpreter','latex','location','Northeast','Fontsize',10);
title('FX rates','interpreter','latex','Fontsize',15);
for j=1:length(curlist)
    temp1 = mean(eval(['ln_' curlist{j} '_us_t((datesperiod))']));
    temp2 = mean(eval(['log(inv_e_' curlist{j} '_t)']));
    subplot(3,3,j+1);
    plot(dates(datesperiod),eval(['-log(inv_e_' curlist{j} '_t)'])+temp2,'LineWidth',2);hold on;
    plot(dates(datesperiod),eval(['ln_' curlist{j} '_us_t((datesperiod))'])-temp1,'LineWidth',2);hold off;
    datetick('x','yy');
    grid on; axis tight;
    label_x('Time (Year)');
    formataxis(gca);
    title([conlist{j} '/USD'],'interpreter','latex','Fontsize',15);
end

figure('Name','FX Devaluation Rates','NumberTitle','off')
subplot(3,3,1);
datesperiod_f=datesperiod(2:end);
datesperiod_l=datesperiod(1:end-1);
temp1=(-log(inv_e_t(datesperiod_f))+log(inv_e_t(datesperiod_l)))*abs_scale/100;
temp2=(ln_eu_us_t(datesperiod_f)-ln_eu_us_t(datesperiod_l))*abs_scale/100;
plot(dates(datesperiod_f),temp1,'LineWidth',2);hold on;
plot(dates(datesperiod_f),temp2,'LineWidth',1,'LineStyle','--');hold off;
datetick('x','yy');
grid on; axis tight;
label_x('Time (Year)');
formataxis(gca);
h=legend('Model','Data');
set(h,'interpreter','latex','location','Northeast','Fontsize',10);
title('FX rates','interpreter','latex','Fontsize',15);
for j=1:length(curlist)
    temp1 = eval(['-log(inv_e_' curlist{j} '_t(datesperiod_f))' ' +log(inv_e_' curlist{j} '_t(datesperiod_l))'])*abs_scale/100;
    temp2 = eval(['ln_' curlist{j} '_us_t(datesperiod_f)' '-ln_' curlist{j} '_us_t(datesperiod_l)'])*abs_scale/100;
    subplot(3,3,j+1);
    plot(dates(datesperiod_f),temp1,'LineWidth',2);hold on;
    plot(dates(datesperiod_f),temp2,'LineWidth',1,'LineStyle','--');hold off;
    datetick('x','yy');
    grid on; axis tight;
    label_x('Time (Year)');
    formataxis(gca);
    title([conlist{j} '/USD corr: ' num2str(corr(temp1,temp2))],'interpreter','latex','Fontsize',15);
end
 
% figure('Name','Sigma Ted','NumberTitle','off')
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(sigma_us_t),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),sigma_us_TED_t(datesperiod),'LineWidth',3); 
% grid on; axis tight;
% datetick('x','yyyy-mm','keeplimits');
% label_x('Time (Year-Month)');
% formataxis(gca);
% orient landscape; 
% h=title('RW - Estimated sequence of $\sigma^*$ (US Ted)');
% set(h,'interpreter','latex','fontsize',20);
% 
% % figure('Name','Sigma BP','NumberTitle','off')
% % plot(dates(datesperiod),ones(1,length(datesperiod))*mean(sigma_us_t),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% % plot(dates(datesperiod),sigma_us_bp_t(datesperiod),'LineStyle',':','LineWidth',2); 
% % % plot(dates(datesperiod),sigma_us_TED_t(datesperiod),'LineWidth',3); 
% % grid on; axis tight;
% % datetick('x','yyyy-mm','keeplimits');
% % label_x('Time (Year-Month)');
% % formataxis(gca);
% % orient landscape; % legend(['average'],'BP (target)','Ted (target)');
% % h=title(['RW - Estimated sequence of $\sigma^* (bp)$']);
% % set(h,'interpreter','latex','fontsize',20);
% 
% % figure('Name','Log Sigma Comparisons','NumberTitle','off')
% % plot(dates(datesperiod),ones(1,length(datesperiod))*log(mean(sigma_us_t)),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% % % plot(dates(datesperiod),log(sigma_us_bp_t(datesperiod)),'LineStyle',':','LineWidth',2); 
% % plot(dates(datesperiod),log(sigma_us_TED_t(datesperiod)),'LineWidth',3); 
% % grid on; axis tight;
% % datetick('x','yyyy-mm','keeplimits');
% % label_x('Time (Year-Month)');
% % formataxis(gca);
% % orient landscape; legend(['average'],'BP (target)','Ted (target)');
% % h=title(['RW - Estimated sequence of log $\sigma^* (bp)$']);
% % set(h,'interpreter','latex','fontsize',20);
% 
% % figure('Name','Sigma EU Comparisons','NumberTitle','off')
% % plot(dates(datesperiod),ones(1,length(datesperiod))*mean(sigma_eu_t),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% % % plot(dates(datesperiod),sigma_eu_t(datesperiod),'LineWidth',3);
% % % plot(dates(datesperiod),sigma_eu_bp_t(datesperiod),'LineWidth',3);
% % plot(dates(datesperiod),sigma_eu_TED_t(datesperiod),'LineWidth',3); 
% % grid on; axis tight;
% % datetick('x','yyyy-mm','keeplimits');
% % label_x('Time (Year-Month)');
% % formataxis(gca);
% % orient landscape; legend([],'BP','TED');
% % set(gcf, 'PaperUnits', 'normalized');
% % set(gcf, 'PaperPosition', [0 0 1 1]);  
% % h=title('RW - Estimated sequence of $\sigma$');
% % set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','Sigma EU Comparisons','NumberTitle','off')
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(sigma_eu_t),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),sigma_eu_TED_t(datesperiod),'LineWidth',3); 
% grid on; axis tight;
% datetick('x','yyyy-mm','keeplimits');
% label_x('Time (Year-Month)');
% formataxis(gca);
% orient landscape; legend([],'BP','TED');
% set(gcf, 'PaperUnits', 'normalized');
% set(gcf, 'PaperPosition', [0 0 1 1]);  
% h=title('RW - Estimated sequence of $\sigma$ (EU Ted)');
% set(h,'interpreter','latex','fontsize',20);
% 
% % figure('Name','Log Sigma EU Comparisons','NumberTitle','off')
% % plot(dates(datesperiod),ones(1,length(datesperiod))*log(mean(sigma_eu_t)),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% % % plot(dates(datesperiod),sigma_eu_t(datesperiod),'LineWidth',3);
% % plot(dates(datesperiod),log(sigma_eu_bp_t(datesperiod)),'LineWidth',3);
% % plot(dates(datesperiod),log(sigma_eu_TED_t(datesperiod)),'LineWidth',3); 
% % grid on; axis tight;
% % datetick('x','yyyy-mm','keeplimits');
% % label_x('Time (Year-Month)');
% % formataxis(gca);
% % orient landscape; legend(['average'],'BP','TED');
% % set(gcf, 'PaperUnits', 'normalized');
% % set(gcf, 'PaperPosition', [0 0 1 1]);  
% % h=title('RW - Estimated sequence of log $\sigma$');
% % set(h,'interpreter','latex','fontsize',20);
% 
% % figure('Name','Log Sigma (bp) Comparisons','NumberTitle','off')
% % plot(dates(datesperiod),ones(1,length(datesperiod))*log(mean(sigma_us_t)),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% % plot(dates(datesperiod),log(sigma_us_bp_t(datesperiod)),'LineStyle',':','LineWidth',2); 
% % plot(dates(datesperiod),log(sigma_eu_bp_t(datesperiod)),'LineWidth',3); 
% % grid on; axis tight;
% % datetick('x','yyyy-mm','keeplimits');
% % label_x('Time (Year-Month)');
% % formataxis(gca);
% % orient landscape; legend(['average'],'US','EU');
% % h=title(['RW - Estimated sequence of $\sigma^*$ vs. $\sigma$ (bp)']);
% % set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','Log Sigma (ted) Comparisons','NumberTitle','off')
% plot(dates(datesperiod),ones(1,length(datesperiod))*log(mean(sigma_us_t)),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),log(sigma_us_TED_t(datesperiod)),'LineStyle',':','LineWidth',2); 
% plot(dates(datesperiod),log(sigma_eu_TED_t(datesperiod)),'LineWidth',3); 
% grid on; axis tight;
% datetick('x','yyyy-mm','keeplimits');
% label_x('Time (Year-Month)');
% formataxis(gca);
% orient landscape; legend('average','US','EU');
% h=title('RW - Estimated sequence of log $\sigma^*$ vs. $\sigma$ (Ted)');
% set(h,'interpreter','latex','fontsize',20);
% 
% % figure
% % plot(dates(datesperiod),sigma_us_TED_flag(datesperiod),'LineWidth',1,'LineStyle',':','Color','y'); hold on;
% % plot(dates(datesperiod),sigma_eu_flag(datesperiod),'LineWidth',1); 
% % grid on; axis tight;
% % datetick('x','yyyy-mm','keeplimits');
% % label_x('Time (Year-Month)');
% % legend('sigma us','sigma eu','rest')
% % formataxis(gca);
% % orient landscape;
% % h=title('Diagnostics');
% % set(h,'interpreter','latex','fontsize',20);
% 
% figure
% plot(dates(datesperiod),sigma_us_TED_flag(datesperiod),'LineWidth',1,'LineStyle','-','Color','b'); hold on;
% plot(dates(datesperiod),c(datesperiod),'LineWidth',1); 
% grid on; axis tight;
% datetick('x','yyyy-mm','keeplimits');
% label_x('Time (Year-Month)');
% legend('sigma us','sigma eu','rest')
% formataxis(gca);
% orient landscape;
% h=title('Diagnostics');
% set(h,'interpreter','latex','fontsize',20);
% 
% % External Fit
% % Risk Premia Outcomes
% figure('Name','Risk Premia','NumberTitle','off') 
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(riskprm_t)*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),riskprm_t(datesperiod)*abs_scale,'LineWidth',3);
% grid on; axis tight;
% datetick('x','yyyy-mm','keeplimits');
% label_x('Time (Year-Month)');
% formataxis(gca);
% orient landscape;
% set(gcf, 'PaperUnits', 'normalized');
% set(gcf, 'PaperPosition', [0 0 1 1]);  
% h=title('RW - Estimated Risk Premium $\xi$');
% set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','Risk Premia vs. Policy Rate Diff','NumberTitle','off') 
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(riskprm_t)*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),riskprm_t(datesperiod)*abs_scale,'LineWidth',3);
% plot(dates(datesperiod),(exp(im_us(datesperiod))-exp(im_eu(datesperiod)))*abs_scale,'LineWidth',2,'LineStyle',':');
% grid on; axis tight; legend([],'model','rate differentials data')
% datetick('x','yyyy-mm','keeplimits');
% label_x('Time (Year-Month)');
% formataxis(gca);
% orient landscape;
% set(gcf, 'PaperUnits', 'normalized');
% set(gcf, 'PaperPosition', [0 0 1 1]);  
% h=title('RW - $\xi$ vs. rate differentials');
% set(h,'interpreter','latex','fontsize',20);
% 
% %% External Fit
% figure('Name','Model vs Data (US BP)','NumberTitle','off') 
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(BP_us_t(datesperiod))*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),BP_us_t(datesperiod)*abs_scale,'LineWidth',3);
% plot(dates(datesperiod),Rb_Rm(datesperiod)*abs_scale,'LineWidth',2,'LineStyle',':');
% grid on; axis tight;
% datetick('x','yyyy-mm','keeplimits');
% label_x('Time (Year-Month)');
% formataxis(gca); legend([],'model','data')
% orient landscape;
% set(gcf, 'PaperUnits', 'normalized');
% set(gcf, 'PaperPosition', [0 0 1 1]);  
% h=title('RW - Model vs Data BP (US)');
% set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','Model vs Data (EU BP)','NumberTitle','off') 
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(BP_eu_t(datesperiod))*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),BP_eu_t(datesperiod)*abs_scale,'LineWidth',3);
% plot(dates(datesperiod),Rb_Rm_eu(datesperiod)*abs_scale,'LineWidth',2,'LineStyle',':');
% grid on; axis tight;
% datetick('x','yyyy-mm','keeplimits');
% label_x('Time (Year-Month)');
% formataxis(gca); legend([],'model','data')
% orient landscape;
% set(gcf, 'PaperUnits', 'normalized');
% set(gcf, 'PaperPosition', [0 0 1 1]);  
% h=title('RW - Model vs Data BP (EU)');
% set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','Model vs Data TED US','NumberTitle','off') 
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(TED_us_t)*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),TED_us_t(datesperiod)*abs_scale,'LineWidth',3);
% plot(dates(datesperiod),(RLibor_us(datesperiod)-im_us(datesperiod))*abs_scale,'LineWidth',3);
% plot(dates(datesperiod),TED_s_us_t(datesperiod)*abs_scale,'LineWidth',3);
% grid on; axis tight;
% datetick('x','yyyy-mm','keeplimits');
% label_x('Time (Year-Month)');
% formataxis(gca);
% orient landscape; legend([],'model','data (diff)','data (dir)')
% set(gcf, 'PaperUnits', 'normalized');
% set(gcf, 'PaperPosition', [0 0 1 1]);  
% h=title('RW - Ted Spread (US)');
% set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','Model vs Data TED EU','NumberTitle','off') 
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(TED_eu_t)*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),TED_s_eu_t(datesperiod)*abs_scale,'LineWidth',3);
% plot(dates(datesperiod),(RLibor_eu(datesperiod)-im_eu(datesperiod))*abs_scale,'LineWidth',3);
% grid on; axis tight;
% datetick('x','yyyy-mm','keeplimits');
% label_x('Time (Year-Month)');
% formataxis(gca);
% orient landscape;
% set(gcf, 'PaperUnits', 'normalized');
% set(gcf, 'PaperPosition', [0 0 1 1]);  
% h=title('RW - Ted Spread (US)');
% set(h,'interpreter','latex','fontsize',20);
% 
% % figure('Name','Model vs CIP','NumberTitle','off') 
% % plot(dates(datesperiod),ones(1,length(datesperiod))*mean(CIP_t(datesperiod))*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% % plot(dates(datesperiod),CIP_t(datesperiod)*abs_scale,'LineWidth',3);
% % plot(dates(datesperiod),CIP_check_t(datesperiod)*abs_scale,'LineWidth',3);
% % plot(dates(datesperiod),cip(datesperiod)*abs_scale,'LineWidth',2,'LineStyle',':');
% % grid on; axis tight;
% % datetick('x','yyyy-mm','keeplimits');
% % label_x('Time (Year-Month)');
% % formataxis(gca); legend('[]','model','model (check)','data')
% % orient landscape;
% % set(gcf, 'PaperUnits', 'normalized');
% % set(gcf, 'PaperPosition', [0 0 1 1]);  
% % h=title('Model vs. Data CIP Deviation');
% % set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','Model vs CIP','NumberTitle','off') 
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(CIP_t(datesperiod))*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),CIP_t(datesperiod)*abs_scale,'LineWidth',3);
% plot(dates(datesperiod),cip(datesperiod)*abs_scale,'LineWidth',2,'LineStyle',':');
% grid on; axis tight;
% datetick('x','yyyy-mm','keeplimits');
% label_x('Time (Year-Month)');
% formataxis(gca); legend([],'model','data')
% orient landscape;
% set(gcf, 'PaperUnits', 'normalized');
% set(gcf, 'PaperPosition', [0 0 1 1]);  
% h=title('Model vs. Data CIP Deviation');
% set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','FX vs. Forward','NumberTitle','off') 
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(CIP_t(datesperiod)),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),exp(-inv_e(datesperiod)),'LineWidth',3);
% plot(dates(datesperiod),f_t(datesperiod),'LineWidth',2,'LineStyle',':');
% grid on; axis tight;
% datetick('x','yyyy-mm','keeplimits');
% label_x('Time (Year-Month)');
% label_y('Level');
% formataxis(gca); legend([],'forward (model)','fx (data)')
% orient landscape;
% set(gcf, 'PaperUnits', 'normalized');
% set(gcf, 'PaperPosition', [0 0 1 1]);  
% h=title('FX vs. Forward');
% set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','Predictable Returns (forward)','NumberTitle','off') 
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(CIP_t(datesperiod))*0,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),(f_t(datesperiod)./exp(-inv_e(datesperiod))-1)*abs_scale,'LineWidth',3);
% plot(dates(datesperiod),(riskprm_t)*abs_scale,'LineWidth',2,'Color','k');
% grid on; axis tight;
% datetick('x','yyyy-mm','keeplimits');
% label_x('Time (Year-Month)');
% label_y('BPS (annual)');
% formataxis(gca); legend([],'forward return','risk premium')
% orient landscape;
% set(gcf, 'PaperUnits', 'normalized');
% set(gcf, 'PaperPosition', [0 0 1 1]);  
% h=title('FX predictable returns');
% set(h,'interpreter','latex','fontsize',20);
% 
% %% Quantity Targets
% % theta_us_t(tt),psi_us_t(tt),Smin_us_t(tt),DW_us_t(tt),FF_us_t(tt)
% figure('Name','Tightness','NumberTitle','off')  
% plot(dates(datesperiod),log(theta_us_t(datesperiod)),'LineWidth',3); hold on;
% plot(dates(datesperiod),log(theta_eu_t(datesperiod)),'LineWidth',3); 
% grid on; axis tight;
% datetick('x','yyyy-mm','keeplimits');
% label_x('Time (Year-Month)');
% formataxis(gca);
% orient landscape;
% set(gcf, 'PaperUnits', 'normalized');
% set(gcf, 'PaperPosition', [0 0 1 1]); 
% legend('US','EU');
% h=title('Market Tightness');
% set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','Funding Changes','NumberTitle','off')  
% plot(dates(datesperiod),(Smin_us_t(datesperiod)),'LineWidth',3); hold on;
% plot(dates(datesperiod),(DW_us_t(datesperiod)),'LineWidth',3); 
% plot(dates(datesperiod),(FF_us_t(datesperiod)),'LineWidth',1); 
% grid on; axis tight;
% datetick('x','yyyy-mm','keeplimits');
% label_x('Time (Year-Month)');
% formataxis(gca);
% orient landscape;
% set(gcf, 'PaperUnits', 'normalized');
% set(gcf, 'PaperPosition', [0 0 1 1]); 
% legend('Funding Need','DW Funded','FF Funded');
% h=title('FX');
% set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','Funding Changes (EU)','NumberTitle','off')  
% plot(dates(datesperiod),(Smin_eu_t(datesperiod)),'LineWidth',3); hold on;
% plot(dates(datesperiod),(DW_eu_t(datesperiod)),'LineWidth',3); 
% plot(dates(datesperiod),(FF_eu_t(datesperiod)),'LineWidth',1);
% grid on; axis tight;
% datetick('x','yyyy-mm','keeplimits');
% label_x('Time (Year-Month)');
% formataxis(gca);
% orient landscape;
% set(gcf, 'PaperUnits', 'normalized');
% set(gcf, 'PaperPosition', [0 0 1 1]); 
% legend('Funding Need','DW Funded','FF Funded');
% h=title('FX');
% set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','Funding (% Dev SS)','NumberTitle','off')  
% plot(dates(datesperiod),(Smin_us_t(datesperiod)/mean(Smin_us_t(datesperiod))),'LineWidth',3); hold on;
% plot(dates(datesperiod),(DW_us_t(datesperiod)/mean(DW_us_t(datesperiod))),'LineWidth',2); 
% plot(dates(datesperiod),(FF_us_t(datesperiod)/mean(FF_us_t(datesperiod))),'LineWidth',1); 
% grid on; axis tight;
% datetick('x','yyyy-mm','keeplimits');
% label_x('Time (Year-Month)');
% formataxis(gca);
% orient landscape;
% set(gcf, 'PaperUnits', 'normalized');
% set(gcf, 'PaperPosition', [0 0 1 1]); 
% legend('Funding Need','DW Funded','FF_us_t');
% h=title('FX');
% set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','Funding EU (% Dev SS)','NumberTitle','off')  
% plot(dates(datesperiod),(Smin_eu_t(datesperiod)/mean(Smin_eu_t(datesperiod))),'LineWidth',3); hold on;
% plot(dates(datesperiod),(DW_eu_t(datesperiod)/mean(DW_eu_t(datesperiod))),'LineWidth',2); 
% plot(dates(datesperiod),(FF_eu_t(datesperiod)/mean(FF_eu_t(datesperiod))),'LineWidth',1); 
% grid on; axis tight;
% datetick('x','yyyy-mm','keeplimits');
% label_x('Time (Year-Month)');
% formataxis(gca);
% orient landscape;
% set(gcf, 'PaperUnits', 'normalized');
% set(gcf, 'PaperPosition', [0 0 1 1]); 
% legend('Funding Need','DW Funded','FF Funded');
% h=title('FX');
% set(h,'interpreter','latex','fontsize',20);
% 
% % Tightness and Interbank Dispersion
% figure('Name','Dispersion (US)','NumberTitle','off') 
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(Chi_D_US(datesperiod))*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),Q_us_t(datesperiod)*abs_scale,'LineWidth',3);
% plot(dates(datesperiod),Chi_D_US(datesperiod)*abs_scale,'LineWidth',2,'LineStyle',':');
% grid on; axis tight;
% datetick('x','yyyy-mm','keeplimits');
% label_x('Time (Year-Month)');
% formataxis(gca); legend([],'model','(Quantiles) data')
% orient landscape;
% set(gcf, 'PaperUnits', 'normalized');
% set(gcf, 'PaperPosition', [0 0 1 1]);  
% h=title('RW - Model vs Data Interbank Dispersion (EU)');
% set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','Dispersion (US)','NumberTitle','off') 
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(Chi_D_US(datesperiod))*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),log(theta_us_t(datesperiod)),'LineWidth',3);
% plot(dates(datesperiod),Chi_D_US(datesperiod)*abs_scale,'LineWidth',2,'LineStyle',':');
% grid on; axis tight;
% datetick('x','yyyy-mm','keeplimits');
% label_x('Time (Year-Month)');
% formataxis(gca); legend([],'model (theta)','data')
% orient landscape;
% set(gcf, 'PaperUnits', 'normalized');
% set(gcf, 'PaperPosition', [0 0 1 1]);  
% h=title('RW - Model vs Data Interbank Dispersion (EU)');
% set(h,'interpreter','latex','fontsize',20);
% 
% % Quantity Variables
% figure('Name','nu eu','NumberTitle','off')
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(nu_t),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),nu_t(datesperiod),'LineWidth',3);
% grid on; axis tight;
% datetick('x','yyyy-mm','keeplimits');
% label_x('Time (Year-Month)');
% formataxis(gca);
% orient landscape;
% h=title('RW - Estimated sequence of $\nu$');
% set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','Theta_d US','NumberTitle','off')
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(Theta_d_us_t),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),Theta_d_us_t(datesperiod),'LineWidth',3);
% grid on; axis tight;
% datetick('x','yyyy-mm','keeplimits');
% label_x('Time (Year-Month)');
% formataxis(gca);
% orient landscape;
% h=title('RW - Estimated sequence of $\Theta_d^*$');
% set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','Theta_d EU','NumberTitle','off')
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(Theta_d_eu_t),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),Theta_d_eu_t(datesperiod),'LineWidth',3);
% grid on; axis tight;
% datetick('x','yyyy-mm','keeplimits');
% label_x('Time (Year-Month)');
% formataxis(gca);
% h=title('RW - Estimated sequence of $\Theta_d$');
% set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','Sigma US (all targets)','NumberTitle','off')  
% plot(dates(datesperiod),0*ones(1,length(datesperiod)),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),log(sigma_us_t(datesperiod)./mean(sigma_us_t(datesperiod))),'LineWidth',3); 
% plot(dates(datesperiod),log(sigma_us_TED_t(datesperiod)./mean(sigma_us_TED_t(datesperiod))),'LineWidth',2); hold on; 
% grid on; axis tight;
% datetick('x','yyyy-mm','keeplimits');
% label_x('Time (Year-Month)');
% formataxis(gca);
% orient landscape;  legend([],'BP','TED');
% h=title('RW - Estimated sequence of $\sigma^*$');
% set(h,'interpreter','latex','fontsize',20);
% 
% % figure('Name','Sigma EU (all targets)','NumberTitle','off')  
% % plot(dates(datesperiod),0*ones(1,length(datesperiod)),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% % plot(dates(datesperiod),log(sigma_eu_t(datesperiod)./mean(sigma_eu_t(datesperiod))),'LineWidth',3); 
% % plot(dates(datesperiod),log(sigma_eu_bp_t(datesperiod)./mean(sigma_eu_bp_t(datesperiod))),'LineWidth',2);
% % plot(dates(datesperiod),log(sigma_eu_TED_t(datesperiod)./mean(sigma_eu_TED_t(datesperiod))),'LineWidth',2); hold on; 
% % grid on; axis tight;
% % datetick('x','yyyy-mm','keeplimits');
% % label_x('Time (Year-Month)');
% % formataxis(gca);
% % orient landscape;  legend('average','CIP','BP','TED');
% % h=title('RW - Estimated sequence of $\sigma^*$');
% % set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','Sigma US (TED target)','NumberTitle','off')  
% plot(dates(datesperiod),0*ones(1,length(datesperiod)),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),log(sigma_us_t(datesperiod)./mean(sigma_us_t(datesperiod))),'LineWidth',3); hold on; 
% plot(dates(datesperiod),log(sigma_eu_t(datesperiod)./mean(sigma_eu_t(datesperiod))),'LineWidth',3); hold on; 
% grid on; axis tight;
% datetick('x','yyyy-mm','keeplimits');
% label_x('Time (Year-Month)');
% formataxis(gca);
% orient landscape; legend([],'US','EU');
% h=title('RW - Estimated sequence of $\sigma^*$');
% set(h,'interpreter','latex','fontsize',20);
% 
% %% Funding Counterfactual
% figure('Name','DW Borrowing (US)','NumberTitle','off') 
% mod_scale  = mean(DW_us_t(datesperiod)./M_us(datesperiod)');
% data_scale = mean(DW_t(datesperiod)./M_us(datesperiod))  ;
% plot(dates(datesperiod),(DW_us_t(datesperiod)./M_us(datesperiod)')/mod_scale,'LineWidth',3); hold on;
% plot(dates(datesperiod),(DW_t(datesperiod)./M_us(datesperiod))/data_scale,'LineWidth',2,'LineStyle',':');
% grid on; axis tight;
% datetick('x','yyyy-mm','keeplimits');
% label_x('Time (Year-Month)');
% formataxis(gca); legend('model','data')
% orient landscape;
% set(gcf, 'PaperUnits', 'normalized');
% set(gcf, 'PaperPosition', [0 0 1 1]);  
% h=title(['RW - Model vs Data DW Borrowing/Reserves']);
% set(h,'interpreter','latex','fontsize',20);
