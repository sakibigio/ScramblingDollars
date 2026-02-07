%% RW_filter_plot_baseline
% Set output folder based on machine
[~, username] = system('whoami');
username = strtrim(username);
if strcmp(username, 'sakibigio')
    foldername = '/Users/sakibigio/Dropbox/Apps/Overleaf/ScramblingDollars_Revision_Restud/quantfigs/';
elseif strcmp(username, 'sakiclaudia')
    foldername = '/Users/sakiclaudia/Library/CloudStorage/Dropbox/Apps/Overleaf/ScramblingDollarsLiquidity_NewVersion/quantfigs/';
else
    error('Unknown user: %s. Please set foldername manually.', username);
end
close all;

% plot baseline results from filtering excercise
figure('Name','Sigma Ted','NumberTitle','off')
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(sigma_us_t),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),sigma_us_TED_t(datesperiod),'LineWidth',3); 
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
orient landscape; 
if printit==1
    exportfig(gcf,[foldername 'F_sigmaus'],'color','cmyk','resolution',1600);
end
h=title('RW - Estimated sequence of $\sigma^*$ (US Ted)');
set(h,'interpreter','latex','fontsize',20);

% figure('Name','Sigma BP','NumberTitle','off')
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(sigma_us_t),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),sigma_us_bp_t(datesperiod),'LineStyle',':','LineWidth',2); 
% % plot(dates(datesperiod),sigma_us_TED_t(datesperiod),'LineWidth',3); 
% grid on; axis tight;
% datetick('x','mmm-yy','keeplimits');
% % label_x('Time (Year-Month)');
% formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
% orient landscape; % legend(['average'],'BP (target)','Ted (target)');
% h=title(['RW - Estimated sequence of $\sigma^* (bp)$']);
% set(h,'interpreter','latex','fontsize',20);

figure('Name','Log Sigma Comparisons','NumberTitle','off')
plot(dates(datesperiod),ones(1,length(datesperiod))*log(mean(sigma_us_t)),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),log(sigma_us_bp_t(datesperiod)),'LineStyle',':','LineWidth',2); 
plot(dates(datesperiod),log(sigma_us_TED_t(datesperiod)),'LineWidth',3); 
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces)); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
orient landscape; legend(['average'],'BP (target)','Ted (target)');
if printit==1
    exportfig(gcf,[foldername 'F_sigmaus_comps'],'color','cmyk','resolution',1600);
end
% h=title(['RW - Estimated sequence of log $\sigma^* (bp)$']);
% set(h,'interpreter','latex','fontsize',20);

% figure('Name','Sigma EU Comparisons','NumberTitle','off')
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(sigma_eu_t),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% % plot(dates(datesperiod),sigma_eu_t(datesperiod),'LineWidth',3);
% % plot(dates(datesperiod),sigma_eu_bp_t(datesperiod),'LineWidth',3);
% plot(dates(datesperiod),sigma_eu_TED_t(datesperiod),'LineWidth',3); 
% grid on; axis tight;
% datetick('x','mmm-yy','keeplimits');
% % label_x('Time (Year-Month)');
% formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
% orient landscape; legend([],'BP','TED');
% set(gcf, 'PaperUnits', 'normalized');
% set(gcf, 'PaperPosition', [0 0 1 1]);  
% h=title('RW - Estimated sequence of $\sigma$');
% set(h,'interpreter','latex','fontsize',20);

figure('Name','Sigma EU Ted','NumberTitle','off')
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(sigma_eu_t),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),sigma_eu_TED_t(datesperiod),'LineWidth',3); 
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
orient landscape; legend([],'BP','TED');
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]);
if printit==1
    exportfig(gcf,[foldername 'F_sigmaeu'],'color','cmyk','resolution',1600);
end
%h=title('RW - Estimated sequence of $\sigma$ (EU Ted)');
% set(h,'interpreter','latex','fontsize',20);

% figure('Name','Log Sigma EU Comparisons','NumberTitle','off')
% plot(dates(datesperiod),ones(1,length(datesperiod))*log(mean(sigma_eu_t)),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% % plot(dates(datesperiod),sigma_eu_t(datesperiod),'LineWidth',3);
% plot(dates(datesperiod),log(sigma_eu_bp_t(datesperiod)),'LineWidth',3);
% plot(dates(datesperiod),log(sigma_eu_TED_t(datesperiod)),'LineWidth',3); 
% grid on; axis tight;
% datetick('x','mmm-yy','keeplimits');
% % label_x('Time (Year-Month)');
% formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
% orient landscape; legend(['average'],'BP','TED');
% set(gcf, 'PaperUnits', 'normalized');
% set(gcf, 'PaperPosition', [0 0 1 1]);  
% h=title('RW - Estimated sequence of log $\sigma$');
% set(h,'interpreter','latex','fontsize',20);

% figure('Name','Log Sigma (bp) Comparisons','NumberTitle','off')
% plot(dates(datesperiod),ones(1,length(datesperiod))*log(mean(sigma_us_t)),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),log(sigma_us_bp_t(datesperiod)),'LineStyle',':','LineWidth',2); 
% plot(dates(datesperiod),log(sigma_eu_bp_t(datesperiod)),'LineWidth',3); 
% grid on; axis tight;
% datetick('x','mmm-yy','keeplimits');
% % label_x('Time (Year-Month)');
% formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
% orient landscape; legend(['average'],'US','EU');
% h=title(['RW - Estimated sequence of $\sigma^*$ vs. $\sigma$ (bp)']);
% set(h,'interpreter','latex','fontsize',20);

figure('Name','Log Sigma (ted) Comparisons','NumberTitle','off')
plot(dates(datesperiod),ones(1,length(datesperiod))*log(mean(sigma_us_t)),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),log(sigma_us_TED_t(datesperiod)),'LineStyle',':','LineWidth',2); 
plot(dates(datesperiod),log(sigma_eu_TED_t(datesperiod)),'LineWidth',3); 
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
orient landscape; legend('average','US','EU');
if printit==1
    exportfig(gcf,[foldername 'F_useu_logcomp'],'color','cmyk','resolution',1600);
end
% h=title('RW - Estimated sequence of log $\sigma^*$ vs. $\sigma$ (Ted)');
% set(h,'interpreter','latex','fontsize',20);

figure('Name','Log Sigma (ted) Comparisons','NumberTitle','off')
plot(dates(datesperiod),ones(1,length(datesperiod))*log(mean(sigma_us_t)),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),(sigma_us_TED_t(datesperiod)),'LineStyle',':','LineWidth',2); 
plot(dates(datesperiod),(sigma_eu_TED_t(datesperiod)),'LineWidth',3); 
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
orient landscape; legend('average','US','EU');
if printit==1
    exportfig(gcf,[foldername 'F_useu_comp'],'color','cmyk','resolution',1600);
end
% h=title('RW - Estimated sequence of log $\sigma^*$ vs. $\sigma$ (Ted)');
% set(h,'interpreter','latex','fontsize',20);

% figure
% plot(dates(datesperiod),sigma_us_TED_flag(datesperiod),'LineWidth',1,'LineStyle',':','Color','y'); hold on;
% plot(dates(datesperiod),sigma_eu_flag(datesperiod),'LineWidth',1); 
% grid on; axis tight;
% datetick('x','mmm-yy','keeplimits');
% % label_x('Time (Year-Month)');
% legend('sigma us','sigma eu','rest')
% formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
% orient landscape;
% h=title('Diagnostics');
% set(h,'interpreter','latex','fontsize',20);

figure
plot(dates(datesperiod),sigma_us_TED_flag(datesperiod),'LineWidth',1,'LineStyle','-','Color','b'); hold on;
plot(dates(datesperiod),sigma_eu_TED_flag(datesperiod),'LineWidth',1); 
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
legend('sigma us','sigma eu','rest')
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
orient landscape;
h=title('Diagnostics');
set(h,'interpreter','latex','fontsize',20);

% External Fit
% Risk Premia Outcomes
figure('Name','Risk Premia','NumberTitle','off') 
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(riskprm_t)*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),riskprm_t(datesperiod)*abs_scale,'LineWidth',3);
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
orient landscape;
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]);  
h=title('RW - Estimated Risk Premium $\xi$');
set(h,'interpreter','latex','fontsize',20);

figure('Name','Risk Premia vs. Policy Rate Diff','NumberTitle','off') 
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(riskprm_t)*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),riskprm_t(datesperiod)*abs_scale,'LineWidth',3);
plot(dates(datesperiod),(exp(im_us(datesperiod))-exp(im_eu(datesperiod)))*abs_scale,'LineWidth',2,'LineStyle',':');
grid on; axis tight; legend([],'model','rate differentials data')
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
orient landscape;
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]);  
h=title('RW - $\xi$ vs. rate differentials');
set(h,'interpreter','latex','fontsize',20);

%% External Fit
figure('Name','Model vs Data (US BP)','NumberTitle','off') 
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(BP_us_t(datesperiod))*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),BP_us_t(datesperiod)*abs_scale,'LineWidth',3);
plot(dates(datesperiod),Rb_Rm(datesperiod)*abs_scale,'LineWidth',2,'LineStyle',':');
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces)); legend([],'model','data')
orient landscape;
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]);  
h=title('RW - Model vs Data BP (US)');
set(h,'interpreter','latex','fontsize',20);

figure('Name','Model vs Data (US BP)','NumberTitle','off') 
plot(dates(datesperiod),zeros(1,length(datesperiod))*mean(BP_us_t(datesperiod))*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),(BP_us_t(datesperiod)-mean(BP_us_t))*abs_scale,'LineWidth',3);
plot(dates(datesperiod),(Rb_Rm(datesperiod)-mean(Rb_Rm(datesperiod)))*abs_scale,'LineWidth',2,'LineStyle',':');
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces)); legend([],'model','data')
orient landscape;
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]);  
h=title('RW - Model vs Data BP (US)');
set(h,'interpreter','latex','fontsize',20);

figure('Name','Model vs Data (EU BP)','NumberTitle','off') 
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(BP_eu_t(datesperiod))*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),BP_eu_t(datesperiod)*abs_scale,'LineWidth',3);
plot(dates(datesperiod),Rb_Rm_eu(datesperiod)*abs_scale,'LineWidth',2,'LineStyle',':');
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces)); legend([],'model','data')
orient landscape;
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]);  
h=title('RW - Model vs Data BP (EU)');
set(h,'interpreter','latex','fontsize',20);

figure('Name','Model vs Data (US BP)','NumberTitle','off') 
plot(dates(datesperiod),zeros(1,length(datesperiod))*mean(BP_eu_t(datesperiod))*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),(BP_eu_t(datesperiod)-mean(BP_eu_t))*abs_scale,'LineWidth',3);
plot(dates(datesperiod),(Rb_Rm_eu(datesperiod)-mean(Rb_Rm_eu(datesperiod)))*abs_scale,'LineWidth',2,'LineStyle',':');
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces)); legend([],'model','data')
orient landscape;
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]);  
h=title('RW - Model vs Data BP (US)');
set(h,'interpreter','latex','fontsize',20);

figure('Name','Model vs Data TED US','NumberTitle','off') 
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(TED_us_t)*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),TED_us_t(datesperiod)*abs_scale,'LineWidth',3);
plot(dates(datesperiod),(RLibor_us(datesperiod)-im_us(datesperiod))*abs_scale,'LineWidth',3);
plot(dates(datesperiod),TED_s_us_t(datesperiod)*abs_scale,'LineWidth',3);
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
orient landscape; legend([],'model','data (diff)','data (dir)')
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]);  
h=title('RW - Ted Spread (US)');
set(h,'interpreter','latex','fontsize',20);

figure('Name','Model vs Data TED EU','NumberTitle','off') 
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(TED_eu_t)*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),TED_s_eu_t(datesperiod)*abs_scale,'LineWidth',3);
plot(dates(datesperiod),(RLibor_eu(datesperiod)-im_eu(datesperiod))*abs_scale,'LineWidth',3);
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
orient landscape;
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]);  
h=title('RW - Ted Spread (US)');
set(h,'interpreter','latex','fontsize',20);

% figure('Name','Model vs CIP','NumberTitle','off') 
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(CIP_t(datesperiod))*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),CIP_t(datesperiod)*abs_scale,'LineWidth',3);
% plot(dates(datesperiod),CIP_check_t(datesperiod)*abs_scale,'LineWidth',3);
% plot(dates(datesperiod),cip(datesperiod)*abs_scale,'LineWidth',2,'LineStyle',':');
% grid on; axis tight;
% datetick('x','mmm-yy','keeplimits');
% % label_x('Time (Year-Month)');
% formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces)); legend('[]','model','model (check)','data')
% orient landscape;
% set(gcf, 'PaperUnits', 'normalized');
% set(gcf, 'PaperPosition', [0 0 1 1]);  
% h=title('Model vs. Data CIP Deviation');
% set(h,'interpreter','latex','fontsize',20);

figure('Name','Model vs CIP','NumberTitle','off') 
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(CIP_t(datesperiod))*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),CIP_t(datesperiod)*abs_scale,'LineWidth',3);
plot(dates(datesperiod),cip(datesperiod)*abs_scale,'LineWidth',2,'LineStyle',':');
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces)); legend([],'model','data')
orient landscape;
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]);  
h=title('Model vs. Data CIP Deviation');
set(h,'interpreter','latex','fontsize',20);

figure('Name','FX vs. Forward','NumberTitle','off') 
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(CIP_t(datesperiod)),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),exp(-inv_e(datesperiod)),'LineWidth',3);
plot(dates(datesperiod),f_t(datesperiod),'LineWidth',2,'LineStyle',':');
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
label_y('Level');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces)); legend([],'forward (model)','fx (data)')
orient landscape;
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]);  
h=title('FX vs. Forward');
set(h,'interpreter','latex','fontsize',20);

figure('Name','Predictable Returns (forward)','NumberTitle','off') 
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(CIP_t(datesperiod))*0,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),(f_t(datesperiod)./exp(-inv_e(datesperiod))-1)*abs_scale,'LineWidth',3);
plot(dates(datesperiod),(riskprm_t)*abs_scale,'LineWidth',2,'Color','k');
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
label_y('BPS (annual)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces)); legend([],'forward return','risk premium')
orient landscape;
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]);  
h=title('FX predictable returns');
set(h,'interpreter','latex','fontsize',20);

%% Quantity Targets
% theta_us_t(tt),psi_us_t(tt),Smin_us_t(tt),DW_us_t(tt),FF_us_t(tt)
figure('Name','Tightness','NumberTitle','off')  
plot(dates(datesperiod),log(theta_us_t(datesperiod)),'LineWidth',3); hold on;
plot(dates(datesperiod),log(theta_eu_t(datesperiod)),'LineWidth',3); 
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
orient landscape;
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]); 
legend('US','EU');
h=title('Market Tightness');
set(h,'interpreter','latex','fontsize',20);

figure('Name','Funding Changes','NumberTitle','off')  
plot(dates(datesperiod),(Smin_us_t(datesperiod)),'LineWidth',3); hold on;
plot(dates(datesperiod),(DW_us_t(datesperiod)),'LineWidth',3); 
plot(dates(datesperiod),(FF_us_t(datesperiod)),'LineWidth',1); 
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
orient landscape;
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]); 
legend('Funding Need','DW Funded','FF Funded');
h=title('FX');
set(h,'interpreter','latex','fontsize',20);

figure('Name','Funding Changes (EU)','NumberTitle','off')  
plot(dates(datesperiod),(Smin_eu_t(datesperiod)),'LineWidth',3); hold on;
plot(dates(datesperiod),(DW_eu_t(datesperiod)),'LineWidth',3); 
plot(dates(datesperiod),(FF_eu_t(datesperiod)),'LineWidth',1);
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
orient landscape;
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]); 
legend('Funding Need','DW Funded','FF Funded');
h=title('FX');
set(h,'interpreter','latex','fontsize',20);

figure('Name','Funding (% Dev SS)','NumberTitle','off')  
plot(dates(datesperiod),(Smin_us_t(datesperiod)/mean(Smin_us_t(datesperiod))),'LineWidth',3); hold on;
plot(dates(datesperiod),(DW_us_t(datesperiod)/mean(DW_us_t(datesperiod))),'LineWidth',2); 
plot(dates(datesperiod),(FF_us_t(datesperiod)/mean(FF_us_t(datesperiod))),'LineWidth',1); 
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
orient landscape;
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]); 
legend('Funding Need','DW Funded','FF_us_t');
h=title('FX');
set(h,'interpreter','latex','fontsize',20);

figure('Name','Funding EU (% Dev SS)','NumberTitle','off')  
plot(dates(datesperiod),(Smin_eu_t(datesperiod)/mean(Smin_eu_t(datesperiod))),'LineWidth',3); hold on;
plot(dates(datesperiod),(DW_eu_t(datesperiod)/mean(DW_eu_t(datesperiod))),'LineWidth',2); 
plot(dates(datesperiod),(FF_eu_t(datesperiod)/mean(FF_eu_t(datesperiod))),'LineWidth',1); 
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
orient landscape;
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]); 
legend('Funding Need','DW Funded','FF Funded');
h=title('FX');
set(h,'interpreter','latex','fontsize',20);

% Tightness and Interbank Dispersion
figure('Name','Dispersion (US)','NumberTitle','off') 
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(Chi_D_US(datesperiod))*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),Q_us_t(datesperiod)*abs_scale,'LineWidth',3);
plot(dates(datesperiod),Chi_D_US(datesperiod)*abs_scale,'LineWidth',2,'LineStyle',':');
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces)); legend([],'model','(Quantiles) data')
orient landscape;
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]);  
h=title('RW - Model vs Data Interbank Dispersion (EU)');
set(h,'interpreter','latex','fontsize',20);

figure('Name','Dispersion (US)','NumberTitle','off') 
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(Chi_D_US(datesperiod))*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),log(theta_us_t(datesperiod)),'LineWidth',3);
plot(dates(datesperiod),Chi_D_US(datesperiod)*abs_scale,'LineWidth',2,'LineStyle',':');
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces)); legend([],'model (theta)','data')
orient landscape;
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]);  
h=title('RW - Model vs Data Interbank Dispersion (EU)');
set(h,'interpreter','latex','fontsize',20);

% Quantity Variables
figure('Name','nu eu','NumberTitle','off')
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(nu_t),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),nu_t(datesperiod),'LineWidth',3);
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
orient landscape;
h=title('RW - Estimated sequence of $\nu$');
set(h,'interpreter','latex','fontsize',20);

figure('Name','Theta_d US','NumberTitle','off')
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(Theta_d_us_t),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),Theta_d_us_t(datesperiod),'LineWidth',3);
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
orient landscape;
h=title('RW - Estimated sequence of $\Theta_d^*$');
set(h,'interpreter','latex','fontsize',20);

figure('Name','Theta_d EU','NumberTitle','off')
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(Theta_d_eu_t),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),Theta_d_eu_t(datesperiod),'LineWidth',3);
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
h=title('RW - Estimated sequence of $\Theta_d$');
set(h,'interpreter','latex','fontsize',20);

figure('Name','d US','NumberTitle','off')
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(d_us_t),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),d_us_t(datesperiod),'LineWidth',3);
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
orient landscape;
h=title('RW - Estimated sequence of $d$ (US)');
set(h,'interpreter','latex','fontsize',20);

figure('Name','d EU','NumberTitle','off')
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(d_eu_t),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),d_eu_t(datesperiod),'LineWidth',3);
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
h=title('RW - Estimated sequence of $d$ (EU)');
set(h,'interpreter','latex','fontsize',20);

figure('Name','nu','NumberTitle','off')
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(nu_t),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),nu_t(datesperiod),'LineWidth',3);
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
h=title('RW - Estimated sequence of $\nu$');
set(h,'interpreter','latex','fontsize',20);

figure('Name','nu','NumberTitle','off')
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(log(nu_t)),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),log(nu_t(datesperiod)),'LineWidth',3);
plot(dates(datesperiod),log(M_eu_t(datesperiod)./M_us_t(datesperiod)),'LineWidth',3);
plot(dates(datesperiod),log(inv_e_t),'LineWidth',3);
plot(dates(datesperiod),log(mu_us_t(datesperiod)./mu_eu_t(datesperiod)),'LineWidth',3);
grid on; axis tight; legend('','nu','M component','FX component','mu component');
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
h=title('RW - Estimated sequence of $\nu$');
set(h,'interpreter','latex','fontsize',20);

figure('Name','log(d) (Comparisons)','NumberTitle','off')
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(log(d_us_t)),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),log(d_us_t(datesperiod)),'LineWidth',3);
plot(dates(datesperiod),log(d_eu_t(datesperiod)),'LineWidth',3);
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
orient landscape; legend([],'US','EU');
h=title('RW - Estimated sequence of $d$');
set(h,'interpreter','latex','fontsize',20);

figure('Name','d (comparisons)','NumberTitle','off')
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(d_us_t),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),d_us_t(datesperiod),'LineWidth',3);
aux=(mu_us_t.*d_us_t)+MBS_us_ss;
plot(dates(datesperiod),aux(datesperiod),'LineWidth',3);
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
orient landscape; legend('','US','Substracted');
h=title('RW - Estimated sequence of $d$');
set(h,'interpreter','latex','fontsize',20);

figure('Name','Sigma US (all targets)','NumberTitle','off')  
plot(dates(datesperiod),0*ones(1,length(datesperiod)),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),log(sigma_us_t(datesperiod)./mean(sigma_us_t(datesperiod))),'LineWidth',3); 
plot(dates(datesperiod),log(sigma_us_TED_t(datesperiod)./mean(sigma_us_TED_t(datesperiod))),'LineWidth',2); hold on; 
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
orient landscape;  legend([],'BP','TED');
h=title('RW - Estimated sequence of $\sigma^*$');
set(h,'interpreter','latex','fontsize',20);

% figure('Name','Sigma EU (all targets)','NumberTitle','off')  
% plot(dates(datesperiod),0*ones(1,length(datesperiod)),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),log(sigma_eu_t(datesperiod)./mean(sigma_eu_t(datesperiod))),'LineWidth',3); 
% plot(dates(datesperiod),log(sigma_eu_bp_t(datesperiod)./mean(sigma_eu_bp_t(datesperiod))),'LineWidth',2);
% plot(dates(datesperiod),log(sigma_eu_TED_t(datesperiod)./mean(sigma_eu_TED_t(datesperiod))),'LineWidth',2); hold on; 
% grid on; axis tight;
% datetick('x','mmm-yy','keeplimits');
% % label_x('Time (Year-Month)');
% formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
% orient landscape;  legend('average','CIP','BP','TED');
% h=title('RW - Estimated sequence of $\sigma^*$');
% set(h,'interpreter','latex','fontsize',20);

figure('Name','Sigma US (TED target)','NumberTitle','off')  
plot(dates(datesperiod),0*ones(1,length(datesperiod)),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),log(sigma_us_t(datesperiod)./mean(sigma_us_t(datesperiod))),'LineWidth',3); hold on; 
plot(dates(datesperiod),log(sigma_eu_t(datesperiod)./mean(sigma_eu_t(datesperiod))),'LineWidth',3); hold on; 
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces));
orient landscape; legend([],'US','EU');
h=title('RW - Estimated sequence of $\sigma^*$');
set(h,'interpreter','latex','fontsize',20);

%% Funding Counterfactual
figure('Name','DW Borrowing (US)','NumberTitle','off') 
mod_scale  = mean(DW_us_t(datesperiod)./M_us(datesperiod)');
data_scale = mean(DW_t(datesperiod)./M_us(datesperiod))  ;
plot(dates(datesperiod),(DW_us_t(datesperiod)./M_us(datesperiod)')/mod_scale,'LineWidth',3); hold on;
plot(dates(datesperiod),(DW_t(datesperiod)./M_us(datesperiod))/data_scale,'LineWidth',2,'LineStyle',':');
grid on; axis tight;
datetick('x','mmm-yy','keeplimits');
% label_x('Time (Year-Month)');
formataxis(gca); ytickformat(gca, sprintf('%%.%df', desiredDecimalPlaces)); legend('model','data')
orient landscape;
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]);  
h=title(['RW - Model vs Data DW Borrowing/Reserves']);
set(h,'interpreter','latex','fontsize',20);
