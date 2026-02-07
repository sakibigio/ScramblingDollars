%% RW_filter_plot_regimes
% saving folder
% saving folder
close all;
[~, username] = system('whoami');
username = strtrim(username);
if strcmp(username, 'sakibigio')
    foldername = '/Users/sakibigio/Dropbox/Apps/Overleaf/ScramblingDollars_Revision_Restud/quantfigs/';
elseif strcmp(username, 'sakiclaudia')
    foldername = '/Users/sakiclaudia/Library/CloudStorage/Dropbox/Apps/Overleaf/ScramblingDollarsLiquidity_NewVersion/quantfigs/';
else
    error('Unknown user: %s', username);
end
printit=0;
eu_color=[0.4 0.1 1.0];
data_color=[0.7 0.1 0.4];
% plot after generating regimes% Import single-state shock
threshold=0.5;
stateprob_tab = readtable('MS_sigma_us_prob.csv', 'TreatAsMissing', 'NA');
sigma_us_stateprob = table2array(stateprob_tab); % single-state shock
sigma_us_stateprob = [sigma_us_stateprob(:,1); 0]; % single-state shock
sigma_us_high_t=sigma_us_stateprob<threshold;

% Change Date Format
if isnumeric(dates)
    dates = datetime(dates, 'ConvertFrom', 'datenum');
elseif ~isdatetime(dates)
    dates = datetime(dates);
end
sigma_us_high2low_t=[(sigma_us_high_t(1:end-1)==1).*(sigma_us_high_t(2:end)==0); 0];
sigma_us_low2high_t=[(sigma_us_high_t(1:end-1)==0).*(sigma_us_high_t(2:end)==1); 0];

% Compute Dates 


% plot baseline results from filtering excercise
figure('Name','High Liquidity Regime','NumberTitle','off')
plot(dates(datesperiod),zeros(1,length(datesperiod))*mean(sigma_us_stateprob),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),sigma_us_stateprob(datesperiod),'LineWidth',3,'Color','r'); 
grid on; axis tight;
xtickformat('MMM-yy');
% label_x('Time (Year-Month)');
formataxis(gca);
orient landscape; 
regime_patches(dates,datesperiod,sigma_us_high_t);
plot(dates(datesperiod),sigma_us_stateprob(datesperiod),'LineWidth',2,'Color','r'); 
if printit==1
   % exportfig(gcf,[foldername 'F_sigmaus_states'],'color','cmyk','resolution',1600);
   exportgraphics(gcf, fullfile(foldername, 'F_sigmaus_states.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end
h=title('Estimated prob. Low Funding Risk State');
set(h,'interpreter','latex','fontsize',20);


% plot baseline results from filtering excercise
figure('Name','Sigma Ted','NumberTitle','off')
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(sigma_us_t),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),sigma_us_TED_t(datesperiod),'LineWidth',3,'LineWidth',3,'Color','r'); 
plot(dates(datesperiod),sigma_eu_TED_t(datesperiod),'LineWidth',3,'LineStyle','-.','Color',eu_color); 
grid on; axis tight;

orient landscape; 
regime_patches(dates,datesperiod,sigma_us_high_t);
plot(dates(datesperiod),sigma_us_TED_t(datesperiod),'LineWidth',3,'Color','r'); 
plot(dates(datesperiod),sigma_eu_TED_t(datesperiod),'LineWidth',3,'LineStyle','-.','Color',eu_color); 
legend('','US','EU','box','off','color','none');
xtickformat('MMM-yy');
% label_x('Time (Year-Month)');
formataxis(gca);
if printit==1
    % exportfig(gcf,[foldername 'F_sigmauseu'],'color','cmyk','resolution',1600);
    exportgraphics(gcf, fullfile(foldername, 'F_sigmauseu.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end
h=title('Estimated sequence of $\sigma^*$ (US Ted)');
set(h,'interpreter','latex','fontsize',20);

figure('Name','TED Comparison','NumberTitle','off') 
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(TED_us_t)*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),TED_us_t(datesperiod)*abs_scale,'LineWidth',3,'Color','r');
plot(dates(datesperiod),TED_eu_t(datesperiod)*abs_scale,'LineWidth',3,'LineStyle','-.','Color',eu_color);
grid on; axis tight;
xtickformat('MMM-yy');
formataxis(gca);
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]);  
regime_patches(dates,datesperiod,sigma_us_high_t);
plot(dates(datesperiod),TED_us_t(datesperiod)*abs_scale,'LineWidth',3,'Color','r');
plot(dates(datesperiod),TED_eu_t(datesperiod)*abs_scale,'LineWidth',3,'LineStyle','-.','Color',eu_color);
if printit==1
    % exportfig(gcf,[foldername 'F_sigmauseu'],'color','cmyk','resolution',1600);
    exportgraphics(gcf, fullfile(foldername, 'F_Tedtargets.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end
h=title('TED Spreads');
set(h,'interpreter','latex','fontsize',20);


%% External Fit
NaNindex=(Rb_Rm==0);
Rb_Rm(NaNindex)=NaN;
figure('Name','US Bond Premium','NumberTitle','off') 
plot(dates(datesperiod),zeros(1,length(datesperiod))*mean(BP_us_t(datesperiod))*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),(BP_us_t(datesperiod)-mean(BP_us_t))*abs_scale,'LineWidth',3);
plot(dates(datesperiod),(Rb_Rm(datesperiod)-mean(Rb_Rm(~isnan(Rb_Rm))))*abs_scale,'LineWidth',3,'LineStyle','-.','Color',data_color);
grid on; axis tight;

% legend([],'model','data')
orient landscape;
regime_patches(dates,datesperiod,sigma_us_high_t);
plot(dates(datesperiod),(BP_us_t(datesperiod)-mean(BP_us_t))*abs_scale,'LineWidth',3,'Color','r');
plot(dates(datesperiod),(Rb_Rm(datesperiod)-mean(Rb_Rm(~isnan(Rb_Rm))))*abs_scale,'LineWidth',3,'LineStyle','-.','Color',data_color);
legend('','Model','Data','box','off','color','none');
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  
xtickformat('MMM-yy');
% label_x('Time (Year-Month)');
formataxis(gca); 
if printit==1
     % exportfig(gcf,[foldername 'F_BPus_fit'],'color','cmyk','resolution',1600);
     exportgraphics(gcf, fullfile(foldername, 'F_BPus_fit.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

figure('Name','Model vs Data (EU BP)','NumberTitle','off') 
NaNindex=(Rb_Rm_eu==0);
Rb_Rm_eu(NaNindex)=NaN;
plot(dates(datesperiod),zeros(1,length(datesperiod))*mean(BP_us_t(datesperiod))*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
plot(dates(datesperiod),(BP_eu_t(datesperiod)-mean(BP_eu_t))*abs_scale,'LineWidth',3);
plot(dates(datesperiod),(Rb_Rm_eu(datesperiod)-mean(Rb_Rm_eu(~isnan(Rb_Rm_eu))))*abs_scale,'LineWidth',3,'LineStyle','-.','Color',data_color);
grid on; axis tight;
% legend([],'model','data')
orient landscape;
regime_patches(dates,datesperiod,sigma_us_high_t);
plot(dates(datesperiod),(BP_eu_t(datesperiod)-mean(BP_eu_t))*abs_scale,'LineWidth',3,'Color','r');
plot(dates(datesperiod),(Rb_Rm_eu(datesperiod)-mean(Rb_Rm_eu(~isnan(Rb_Rm_eu))))*abs_scale,'LineWidth',3,'LineStyle','-.','Color',data_color);
xtickformat('MMM-yy');
% label_x('Time (Year-Month)');
formataxis(gca); 
legend('','Model','Data','box','off','color','none');
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  
if printit==1
     % exportfig(gcf,[foldername 'F_BPeu_fit'],'color','cmyk','resolution',1600);
     exportgraphics(gcf, fullfile(foldername, 'F_BPeu_fit.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

figure('Name','US TED Spread','NumberTitle','off') 
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(TED_us_t)*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),TED_us_t(datesperiod)*abs_scale,'LineWidth',3);
%plot(dates(datesperiod),(RLibor_us(datesperiod)-im_us(datesperiod))*abs_scale,'LineWidth',3);
plot(dates(datesperiod),TED_s_us_t(datesperiod)*abs_scale,'LineWidth',3);
grid on; axis tight;
xtickformat('MMM-yy');
% label_x('Time (Year-Month)');
formataxis(gca);
%orient landscape; legend([],'model','data (diff)','data (dir)')
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]);  
h=title('US TED Spread');
set(h,'interpreter','latex','fontsize',20);
regime_patches(dates,datesperiod,sigma_us_high_t);
plot(dates(datesperiod),TED_s_us_t(datesperiod)*abs_scale,'LineWidth',3);

% figure('Name','Model vs CIP','NumberTitle','off') 
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(CIP_t(datesperiod))*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),CIP_t(datesperiod)*abs_scale,'LineWidth',3);
% plot(dates(datesperiod),CIP_check_t(datesperiod)*abs_scale,'LineWidth',3);
% plot(dates(datesperiod),cip(datesperiod)*abs_scale,'LineWidth',2,'LineStyle',':');
% grid on; axis tight;
% xtickformat('MMM-yy');
% % label_x('Time (Year-Month)');
% formataxis(gca); legend('[]','model','model (check)','data')
% orient landscape;
% set(gcf, 'PaperUnits', 'normalized');
% set(gcf, 'PaperPosition', [0 0 1 1]);  
% h=title('Model vs. Data CIP Deviation');
% set(h,'interpreter','latex','fontsize',20);

figure('Name','CIP deviation (Bond basis)','NumberTitle','off') 
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(CIP_t(datesperiod))*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),CIP_t(datesperiod)*abs_scale,'LineWidth',3);
plot(dates(datesperiod),cip(datesperiod)*abs_scale,'LineWidth',3,'LineStyle','-');
grid on; axis tight;
xtickformat('MMM-yy');
% label_x('Time (Year-Month)');
formataxis(gca);
%legend([],'model','data')
orient landscape;
% set(gcf, 'PaperUnits', 'normalized');
% set(gcf, 'PaperPosition', [0 0 1 1]);  
% h=title('Data CIP (bond) Deviation');
% set(h,'interpreter','latex','fontsize',20);
regime_patches(dates,datesperiod,sigma_us_high_t);
plot(dates(datesperiod),cip(datesperiod)*abs_scale,'LineWidth',3,'LineStyle','-','Color',data_color);
if printit==1
     exportfig(gcf,[foldername 'F_CIP'],'color','cmyk','resolution',1600);
end

%% Quantity Targets
% theta_us_t(tt),psi_us_t(tt),Smin_us_t(tt),DW_us_t(tt),FF_us_t(tt)
figure('Name','Tightness','NumberTitle','off')  
plot(dates(datesperiod),log(theta_us_t(datesperiod)),'LineWidth',3); hold on;
plot(dates(datesperiod),log(theta_eu_t(datesperiod)),'LineWidth',3); 
grid on; axis tight;
xtickformat('MMM-yy');
% label_x('Time (Year-Month)');
formataxis(gca);
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
xtickformat('MMM-yy');
% label_x('Time (Year-Month)');
formataxis(gca);
orient landscape;
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]); 
legend('Funding Need','DW Funded','FF Funded');
h=title('FX');
set(h,'interpreter','latex','fontsize',20);

% figure('Name','Funding Changes (EU)','NumberTitle','off')  
% plot(dates(datesperiod),(Smin_eu_t(datesperiod)),'LineWidth',3); hold on;
% plot(dates(datesperiod),(DW_eu_t(datesperiod)),'LineWidth',3); 
% plot(dates(datesperiod),(FF_eu_t(datesperiod)),'LineWidth',1);
% grid on; axis tight;
% xtickformat('MMM-yy');
% % label_x('Time (Year-Month)');
% formataxis(gca);
% orient landscape;
% set(gcf, 'PaperUnits', 'normalized');
% set(gcf, 'PaperPosition', [0 0 1 1]); 
% legend('Funding Need','DW Funded','FF Funded');
% h=title('FX');
% set(h,'interpreter','latex','fontsize',20);

figure('Name','Funding (% Dev SS)','NumberTitle','off')  
plot(dates(datesperiod),(Smin_us_t(datesperiod)/mean(Smin_us_t(datesperiod))),'LineWidth',3); hold on;
plot(dates(datesperiod),(DW_us_t(datesperiod)/mean(DW_us_t(datesperiod))),'LineWidth',2); 
plot(dates(datesperiod),(FF_us_t(datesperiod)/mean(FF_us_t(datesperiod))),'LineWidth',1); 
grid on; axis tight;
xtickformat('MMM-yy');
% label_x('Time (Year-Month)');
formataxis(gca);
orient landscape;
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]); 
h=legend('Funding Need','DW Funded','FF_us_t');
set(h,'interpreter','latex','fontsize',20);

% figure('Name','Funding EU (% Dev SS)','NumberTitle','off')  
% plot(dates(datesperiod),(Smin_eu_t(datesperiod)/mean(Smin_eu_t(datesperiod))),'LineWidth',3); hold on;
% plot(dates(datesperiod),(DW_eu_t(datesperiod)/mean(DW_eu_t(datesperiod))),'LineWidth',2); 
% plot(dates(datesperiod),(FF_eu_t(datesperiod)/mean(FF_eu_t(datesperiod))),'LineWidth',1); 
% grid on; axis tight;
% xtickformat('MMM-yy');
% % label_x('Time (Year-Month)');
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
% xtickformat('MMM-yy');
% % label_x('Time (Year-Month)');
% formataxis(gca); legend([],'model','(Quantiles) data')
% orient landscape;
% set(gcf, 'PaperUnits', 'normalized');
% set(gcf, 'PaperPosition', [0 0 1 1]);  
% h=title('RW - Model vs Data Interbank Dispersion (EU)');
% set(h,'interpreter','latex','fontsize',20);

figure('Name','Dispersion (US)','NumberTitle','off') 
plot(dates(datesperiod),ones(1,length(datesperiod))*mean(Chi_D_US(datesperiod))*abs_scale,'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),log(theta_us_t(datesperiod)),'LineWidth',3);
plot(dates(datesperiod),Chi_D_US(datesperiod)*abs_scale,'LineWidth',3,'LineStyle','-');
grid on; axis tight;
xtickformat('MMM-yy');
% label_x('Time (Year-Month)');
formataxis(gca); 
%legend([],'model (theta)','data')

% h=title('US Interbank Dispersion');
% set(h,'interpreter','latex','fontsize',20);
regime_patches(dates,datesperiod,sigma_us_high_t)
plot(dates(datesperiod),Chi_D_US(datesperiod)*abs_scale,'LineWidth',2,'LineStyle','-','Color',data_color);
xtickformat('MMM-yy');
% label_x('Time (Year-Month)');
formataxis(gca); 
if printit==1
    orient landscape;
    set(gcf, 'PaperUnits', 'normalized');
    set(gcf, 'PaperPosition', [0 0 1 1]);  
    orient landscape;
     % exportfig(gcf,[foldername 'F_IBdispersion'],'color','cmyk','resolution',1600);
    exportgraphics(gcf, fullfile(foldername, 'F_IBdispersion.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

% % Quantity Variables
% figure('Name','nu eu','NumberTitle','off')
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(nu_t),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),nu_t(datesperiod),'LineWidth',3);
% grid on; axis tight;
% xtickformat('MMM-yy');
% % label_x('Time (Year-Month)');
% formataxis(gca);
% orient landscape;
% h=title('RW - Estimated sequence of $\nu$');
% set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','Theta_d US','NumberTitle','off')
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(Theta_d_us_t),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),Theta_d_us_t(datesperiod),'LineWidth',3);
% grid on; axis tight;
% xtickformat('MMM-yy');
% % label_x('Time (Year-Month)');
% formataxis(gca);
% orient landscape;
% h=title('RW - Estimated sequence of $\Theta_d^*$');
% set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','Theta_d EU','NumberTitle','off')
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(Theta_d_eu_t),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),Theta_d_eu_t(datesperiod),'LineWidth',3);
% grid on; axis tight;
% xtickformat('MMM-yy');
% % label_x('Time (Year-Month)');
% formataxis(gca);
% h=title('RW - Estimated sequence of $\Theta_d$');
% set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','d US','NumberTitle','off')
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(d_us_t),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),d_us_t(datesperiod),'LineWidth',3);
% grid on; axis tight;
% xtickformat('MMM-yy');
% % label_x('Time (Year-Month)');
% formataxis(gca);
% orient landscape;
% h=title('RW - Estimated sequence of $d$ (US)');
% set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','d EU','NumberTitle','off')
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(d_eu_t),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),d_eu_t(datesperiod),'LineWidth',3);
% grid on; axis tight;
% xtickformat('MMM-yy');
% % label_x('Time (Year-Month)');
% formataxis(gca);
% h=title('RW - Estimated sequence of $d$ (EU)');
% set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','nu','NumberTitle','off')
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(nu_t),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),nu_t(datesperiod),'LineWidth',3);
% grid on; axis tight;
% xtickformat('MMM-yy');
% % label_x('Time (Year-Month)');
% formataxis(gca);
% h=title('RW - Estimated sequence of $\nu$');
% set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','nu','NumberTitle','off')
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(log(nu_t)),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),log(nu_t(datesperiod)),'LineWidth',3);
% plot(dates(datesperiod),log(M_eu_t(datesperiod)./M_us_t(datesperiod)),'LineWidth',3);
% plot(dates(datesperiod),log(inv_e_t),'LineWidth',3);
% plot(dates(datesperiod),log(mu_us_t(datesperiod)./mu_eu_t(datesperiod)),'LineWidth',3);
% grid on; axis tight; legend('','nu','M component','FX component','mu component');
% xtickformat('MMM-yy');
% % label_x('Time (Year-Month)');
% formataxis(gca);
% h=title('RW - Estimated sequence of $\nu$');
% set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','log(d) (Comparisons)','NumberTitle','off')
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(log(d_us_t)),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),log(d_us_t(datesperiod)),'LineWidth',3);
% plot(dates(datesperiod),log(d_eu_t(datesperiod)),'LineWidth',3);
% grid on; axis tight;
% xtickformat('MMM-yy');
% % label_x('Time (Year-Month)');
% formataxis(gca);
% orient landscape; legend([],'US','EU');
% h=title('RW - Estimated sequence of $d$');
% set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','d (comparisons)','NumberTitle','off')
% plot(dates(datesperiod),ones(1,length(datesperiod))*mean(d_us_t),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),d_us_t(datesperiod),'LineWidth',3);
% aux=(mu_us_t.*d_us_t)+MBS_us_ss;
% plot(dates(datesperiod),aux(datesperiod),'LineWidth',3);
% grid on; axis tight;
% xtickformat('MMM-yy');
% % label_x('Time (Year-Month)');
% formataxis(gca);
% orient landscape; legend('','US','Substracted');
% h=title('RW - Estimated sequence of $d$');
% set(h,'interpreter','latex','fontsize',20);
% 
% figure('Name','Sigma US (all targets)','NumberTitle','off')  
% plot(dates(datesperiod),0*ones(1,length(datesperiod)),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),log(sigma_us_t(datesperiod)./mean(sigma_us_t(datesperiod))),'LineWidth',3); 
% plot(dates(datesperiod),log(sigma_us_TED_t(datesperiod)./mean(sigma_us_TED_t(datesperiod))),'LineWidth',2); hold on; 
% grid on; axis tight;
% xtickformat('MMM-yy');
% % label_x('Time (Year-Month)');
% formataxis(gca);
% orient landscape;  legend([],'BP','TED');
% h=title('RW - Estimated sequence of $\sigma^*$');
% set(h,'interpreter','latex','fontsize',20);

% figure('Name','Sigma EU (all targets)','NumberTitle','off')  
% plot(dates(datesperiod),0*ones(1,length(datesperiod)),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),log(sigma_eu_t(datesperiod)./mean(sigma_eu_t(datesperiod))),'LineWidth',3); 
% plot(dates(datesperiod),log(sigma_eu_bp_t(datesperiod)./mean(sigma_eu_bp_t(datesperiod))),'LineWidth',2);
% plot(dates(datesperiod),log(sigma_eu_TED_t(datesperiod)./mean(sigma_eu_TED_t(datesperiod))),'LineWidth',2); hold on; 
% grid on; axis tight;
% xtickformat('MMM-yy');
% % label_x('Time (Year-Month)');
% formataxis(gca);
% orient landscape;  legend('average','CIP','BP','TED');
% h=title('RW - Estimated sequence of $\sigma^*$');
% set(h,'interpreter','latex','fontsize',20);

% figure('Name','Sigma US (TED target)','NumberTitle','off')  
% plot(dates(datesperiod),0*ones(1,length(datesperiod)),'LineWidth',2,'LineStyle','--','Color','k'); hold on;
% plot(dates(datesperiod),log(sigma_us_t(datesperiod)./mean(sigma_us_t(datesperiod))),'LineWidth',3); hold on; 
% plot(dates(datesperiod),log(sigma_eu_t(datesperiod)./mean(sigma_eu_t(datesperiod))),'LineWidth',3); hold on; 
% grid on; axis tight;
% xtickformat('MMM-yy');
% % label_x('Time (Year-Month)');
% formataxis(gca);
% orient landscape; legend([],'US','EU');
% h=title('RW - Estimated sequence of $\sigma^*$');
% set(h,'interpreter','latex','fontsize',20);

%% Funding Counterfactual
figure('Name','DW Borrowing (US)','NumberTitle','off') 
mod_scale  = mean(DW_us_t(datesperiod)./M_us(datesperiod)');
data_scale = mean(DW_t(datesperiod)./M_us(datesperiod))  ;
plot(dates(datesperiod),(DW_us_t(datesperiod)./M_us(datesperiod)')/mod_scale,'LineWidth',3,'Color','r');  hold on;
plot(dates(datesperiod),(DW_t(datesperiod)./M_us(datesperiod))/data_scale,'LineWidth',3,'LineStyle','-.','Color',data_color);
grid on; axis tight;
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]);  
orient landscape;
regime_patches(dates,datesperiod,sigma_us_high_t)
plot(dates(datesperiod),(DW_us_t(datesperiod)./M_us(datesperiod)')/mod_scale,'LineWidth',3,'Color','r');  hold on;
plot(dates(datesperiod),(DW_t(datesperiod)./M_us(datesperiod))/data_scale,'LineWidth',3,'LineStyle','-.','Color',data_color);
xtickformat('MMM-yy');
% label_x('Time (Year-Month)');
formataxis(gca); 
legend('Model','Data','box','off','color','none');
set(h,'interpreter','latex','fontsize',20);
if printit==1
    set(gcf, 'PaperUnits', 'normalized');
    set(gcf, 'PaperPosition', [0 0 1 1]);  
    orient landscape;
    %  exportfig(gcf,[foldername 'F_DWus_fit'],'color','cmyk','resolution',1600);
     exportgraphics(gcf, fullfile(foldername, 'F_DWus_fit.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end
%set(gcf, 'PaperUnits', 'normalized');
%set(gcf, 'PaperPosition', [0 0 1 1]);  
%h=title('DW Borrowing/Reserves');
%set(h,'interpreter','latex','fontsize',20);


%% Cross-Country CIP deviations
desiredNumXTicks = 3; % Change this to the number of ticks you want

% CIP Deviations
figure('Name','CIP Deviations','NumberTitle','off') 
subplot(3,3,1); 
% plot(dates(datesperiod),CIP_t(datesperiod)*abs_scale,'LineWidth',2);hold on;
plot(dates(datesperiod),CIP_s_eu_t(datesperiod)*abs_scale,'LineWidth',2);hold on;
xtickformat('MM-yy');
grid on; axis tight;
regime_patches(dates,datesperiod,sigma_us_high_t)
plot(dates(datesperiod),CIP_s_eu_t(datesperiod)*abs_scale,'LineWidth',2,'Color',data_color);hold on;
% label_x('Time (Year)');
ax = gca;
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
% formataxis(gca);
% h=legend('Model','Data');
% set(h,'interpreter','latex','location','Northeast','Fontsize',10);
title('EU','interpreter','latex','Fontsize',15);
for j=1:length(curlist)
    aux = eval(['CIP_' curlist{j} '_t(datesperiod)']);
    aux_d = eval(['CIP_s_' curlist{j} '_t(datesperiod)']);
    subplot(3,3,j+1);
    plot(dates(datesperiod),aux_d(datesperiod)*abs_scale,'LineWidth',2); hold on;
    regime_patches(dates,datesperiod,sigma_us_high_t);
    plot(dates(datesperiod),aux_d(datesperiod)*abs_scale,'LineWidth',2,'Color',data_color);hold on;
    xtickformat('MM-yy');
    grid on; axis tight;
    % label_x('Time (Year)');
    ax = gca;
    xLimits = get(ax, 'XLim');
    customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
    set(ax, 'XTick', customXTicks);
  %  formataxis(gca);
    title([conlist{j}],'interpreter','latex','Fontsize',15);
end
if printit==1
    set(gcf, 'PaperUnits', 'normalized');
    orient landscape;
    set(gcf, 'PaperPosition', [0 0 1 1]);  
    % exportfig(gcf,[foldername 'F_CIP_all'],'color','cmyk','resolution',1600);
    exportgraphics(gcf, fullfile(foldername, 'F_CIP_all.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

figure('Name','FX Devaluation Rates','NumberTitle','off')
subplot(3,3,1);
datesperiod_f=datesperiod(2:end);
datesperiod_l=datesperiod(1:end-1);
temp1=(-log(inv_e_t(datesperiod_f))+log(inv_e_t(datesperiod_l)))*abs_scale/100;
temp2=(ln_eu_us_t(datesperiod_f)-ln_eu_us_t(datesperiod_l))*abs_scale/100;
% plot(dates(datesperiod_f),temp1,'LineWidth',2);hold on;
plot(dates(datesperiod_f),temp2,'LineWidth',1,'LineStyle','-','Color',data_color); hold on;
regime_patches(dates,datesperiod,sigma_us_high_t);
plot(dates(datesperiod_f),temp2,'LineWidth',1,'LineStyle','-','Color',data_color); hold on;
xtickformat('MM-yy');
ax = gca;
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
% formataxis(gca);
grid on; axis tight;
title('EU/USD','interpreter','latex','Fontsize',15);
% label_x('Time (Year)');
% formataxis(gca);
% h=legend('Model','Data');
% set(h,'interpreter','latex','location','Northeast','Fontsize',10);
% title('FX rates','interpreter','latex','Fontsize',15);
for j=1:length(curlist)
    temp1 = eval(['-log(inv_e_' curlist{j} '_t(datesperiod_f))' ' +log(inv_e_' curlist{j} '_t(datesperiod_l))'])*abs_scale/100;
    temp2 = eval(['ln_' curlist{j} '_us_t(datesperiod_f)' '-ln_' curlist{j} '_us_t(datesperiod_l)'])*abs_scale/100;
    subplot(3,3,j+1);
    % plot(dates(datesperiod_f),temp1,'LineWidth',2);hold on;
    plot(dates(datesperiod_f),temp2,'LineWidth',1,'LineStyle','-','Color',data_color); hold on;
    regime_patches(dates,datesperiod,sigma_us_high_t); hold on;
    plot(dates(datesperiod_f),temp2,'LineWidth',1,'LineStyle','-','Color',data_color); hold on;
    xtickformat('MM-yy');
    grid on; axis tight;
    % label_x('Time (Year)');
    ax = gca;
    xLimits = get(ax, 'XLim');
    customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
    set(ax, 'XTick', customXTicks);
    % formataxis(gca);
    title([conlist{j} '/USD' ],'interpreter','latex','Fontsize',15);
end
if printit==1
    set(gcf, 'PaperUnits', 'normalized');
    orient landscape;
    set(gcf, 'PaperPosition', [0 0 1 1]);      
   %  exportfig(gcf,[foldername 'F_devaluation_all'],'color','cmyk','resolution',1600);
    exportgraphics(gcf, fullfile(foldername, 'F_devaluation_all.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

%% Report Regime Moments
% Indexes
index_r2=sigma_us_high_t;
index_r1=~index_r2;
dev_t=((inv_e_t(2:end)./inv_e_t(1:end-1)).^(-1)-1)*abs_scale;

% Average Regime Change
dev_low2high=mean(dev_t(sigma_us_low2high_t==1));
dev_high2low=mean(dev_t(sigma_us_high2low_t==1));

% Unconditional Moments
% means:
E_e_data=mean(inv_e_t);
E_bp_data=mean(BP_us_t)*abs_scale;
E_cip_data=mean(CIP_s_eu_t)*abs_scale;
E_dev_data=mean(dev_t);

% standard deviations:
std_e_data=std(inv_e_t);
std_bp_data=std(BP_us_t)*abs_scale;
std_cip_data=std(CIP_s_eu_t)*abs_scale;
std_dev_data=std(dev_t);

% autocorrelation:
aux=autocorr(inv_e_t);
rho_e_data=aux(2);
aux=autocorr(BP_us_t);
rho_bp_data=aux(2);
aux=autocorr(CIP_s_eu_t);
rho_cip_data=aux(2);
aux=autocorr(dev_t);
rho_dev_data=aux(2);

% Conditional Moments
% means:
E_e_data_r1=mean(inv_e_t(index_r1));
E_bp_data_r1=mean(BP_us_t(index_r1))*abs_scale;
E_cip_data_r1=mean(CIP_s_eu_t(index_r1))*abs_scale;
E_dev_data_r1=mean(dev_t(index_r1(1:end-1)));

% standard deviations:
std_e_data_r1=std(inv_e_t(index_r1));
std_bp_data_r1=std(BP_us_t(index_r1))*abs_scale;
std_cip_data_r1=std(CIP_s_eu_t(index_r1))*abs_scale;
std_dev_data_r1=std(dev_t(index_r1(1:end-1)));

% autocorrelation:
aux=autocorr(inv_e_t(index_r1));
rho_e_data_r1=aux(2);
aux=autocorr(BP_us_t(index_r1));
rho_bp_data_r1=aux(2);
aux=autocorr(CIP_s_eu_t(index_r1));
rho_cip_data_r1=aux(2);
aux=autocorr(dev_t(index_r1(1:end-1)));
rho_dev_data_r1=aux(2);

% Conditional Moments
% means:
E_e_data_r2=mean(inv_e_t(index_r2));
E_bp_data_r2=mean(BP_us_t(index_r2))*abs_scale;
E_cip_data_r2=mean(CIP_s_eu_t(index_r2))*abs_scale;
E_dev_data_r2=mean(dev_t(index_r2(1:end-1)));

% standard deviations:
std_e_data_r2=std(inv_e_t(index_r2));
std_bp_data_r2=std(BP_us_t(index_r2))*abs_scale;
std_cip_data_r2=std(CIP_s_eu_t(index_r2))*abs_scale;
std_dev_data_r2=std(dev_t(index_r2(1:end-1)));

% autocorrelation:
aux=autocorr(inv_e_t(index_r2));
rho_e_data_r2=aux(2);
aux=autocorr(BP_us_t(index_r2));
rho_bp_data_r2=aux(2);
aux=autocorr(CIP_s_eu_t(index_r2));
rho_cip_data_r2=aux(2);
aux=autocorr(dev_t(index_r2(1:end-1)));
rho_dev_data_r2=aux(2);

% Replace the following with your actual model moments
filename = fullfile(foldername, 'Data_Moments.tex');
fid = fopen(filename, 'wt');
data_moments = {'FX', 1, rho_e_data, std_e_data/E_e_data, (E_e_data_r2./E_e_data_r1-1)*10000, rho_e_data_r2, rho_e_data_r1, (std_e_data_r2./std_e_data_r1);
                '$\Delta$ FX', E_dev_data, rho_dev_data, std_dev_data, E_dev_data_r2-E_dev_data_r1, rho_dev_data_r2, rho_dev_data_r1, (std_dev_data_r2./std_dev_data_r1);
                 'BP', E_bp_data, rho_bp_data, std_bp_data, E_bp_data_r2-E_bp_data_r1, rho_bp_data_r2, rho_bp_data_r1, (std_bp_data_r2./std_bp_data_r1);                
                 'CIP', E_cip_data, rho_cip_data, std_cip_data, E_cip_data_r2-E_cip_data_r1, rho_cip_data_r2, rho_cip_data_r1, (std_cip_data_r2./std_cip_data_r1)};

% fprintf(fid, '\\bottomrule\n');
for ii = 1:size(data_moments, 1)
    fprintf(fid, '%s & %.1f & %.2f & %.1f & %.1f &  \\{ %.2f, %.2f \\} & %.1f \\\\ \n', data_moments{ii, :});
end

filename = fullfile(foldername, 'Data_CIP_Moments.tex');
fid = fopen(filename, 'wt');
data_moments = {...
    %'FX', 1, rho_e_data, std_e_data/E_e_data, (E_e_data_r2./E_e_data_r1-1)*10000, rho_e_data_r2, rho_e_data_r1, (std_e_data_r2./std_e_data_r1);
              %  '$\Delta$ FX', E_dev_data, rho_dev_data, std_dev_data, E_dev_data_r2-E_dev_data_r1, rho_dev_data_r2, rho_dev_data_r1, (std_dev_data_r2./std_dev_data_r1);
              %   'BP', E_bp_data, rho_bp_data, std_bp_data, E_bp_data_r2-E_bp_data_r1, rho_bp_data_r2, rho_bp_data_r1, (std_bp_data_r2./std_bp_data_r1);                
                 'CIP (data)', E_cip_data, rho_cip_data, std_cip_data, E_cip_data_r2-E_cip_data_r1, rho_cip_data_r2, rho_cip_data_r1, (std_cip_data_r2./std_cip_data_r1)};

% fprintf(fid, '\\bottomrule\n');
for ii = 1:size(data_moments, 1)
    fprintf(fid, '%s & %.1f & %.2f & %.1f & %.1f &  \\{ %.2f, %.2f \\} & %.1f \\\\ \n', data_moments{ii, :});
end

% Close the table
% fprintf(fid, '\\bottomrule \n');
fclose(fid);



%% Functions
% Define shading for the state = 1 regions
% Loop through the dates and create patches for periods where sigma_us_high_t == 1
function out=regime_patches(dates,datesperiod,sigma_us_high_t)
    % Get the y-axis limits for shading
    y_limits = ylim; % [y_min, y_max]
    for i = 1:length(datesperiod) - 1
        if sigma_us_high_t(i) == 1
            % Define x-coordinates for the patch (spanning two consecutive dates)
            x_patch = [dates(datesperiod(i)), dates(datesperiod(i + 1)), ...
                       dates(datesperiod(i + 1)), dates(datesperiod(i))];
            % Define y-coordinates (spanning the full y-axis range)
            y_patch = [y_limits(1), y_limits(1), y_limits(2), y_limits(2)];
            % Create the shaded patch with semi-transparent color
            patch(x_patch, y_patch, [0.9, 0.9, 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
        end
    end
    out=[];
end

