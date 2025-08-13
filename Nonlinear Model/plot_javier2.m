close all

 
%% Final Plots for paper
% 1) Combine dollar and euro liquidity ratio in one plot [x]
% 2) Combine Rm,Rm* in one plot
% 3) Exchange rate
% 4) BP.
% 5) CIP
% 6) Dollar and Euro deposit rates in one plot

% policy functions
cc = cc+1; 
x_vec = sigma_us_vec(:)/mean(sigma_us_vec);
x_lab = '$\sigma_{us}$';
% x_vec = (im_us_vec.^(freq)-1)*100;
% x_lab = '$i^{us}_m(\%)$';

 num_Tick=3

% desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% % Format the x-axis tick labels
% xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));

foldername='C:\Users\Javier\Dropbox\Apps\Overleaf\ScramblingDollarsLiquidity\newfigs';
 
xticks_def= [0 2 4]

xlim_def=[0 5]

FSize=56
 FS_plot  = 22;
set(0,'defaulttextInterpreter','latex')  
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultLineLineWidth',2,'defaultFigureColor','w')
 
 foldername='';

 col_eu =rgb('red')
% 
%  col_mac   = [ 0     0     1];
col_eu   = [0.70 0.00 0.00];

figure; % Liquidity Ratios
plot(x_vec,mu_us_vec,'b', 'linewidth',2.5); hold on;
plot(x_vec,mu_eu_vec,'color',col_eu, 'linewidth',2.5,'LineStyle','--');
xlim([x_vec(1) x_vec(end)]);
xlab=xlabel(x_lab,'Fontsize',FSize)
% ylabel('$liquidity$','interpreter','latex','Fontsize',FSize);
% %formataxis(gca);
h=legend('$\mu^{us}$','$\mu^{eu}$','color','none','box','off','location','southwest','FontSize',20);
set(h,...
    'Position',[0.710714285711525 0.603968236654524 0.15685206625945 0.163095233553932],...
    'FontSize',20,...
    'Color','none');
xticks(xticks_def)
set(gca,'XLim',xlim_def)
 
% set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', FSize);
% title('Liquidity','interpreter','latex','Fontsize',15);
ax = gca;
% Specify the number of desired ticks on the x-axis
% desiredNumXTicks = num_Tick; % Change this to the number of ticks you want
% xLimits = get(ax, 'XLim');
% customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
% set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);

% desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% % Format the x-axis tick labels
% xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));

 
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
  set(gca,'Fontsize',FS_plot)
% set(gcf, 'PaperUnits', 'normalized');
% orient landscape;
% set(gcf, 'PaperPosition', [0 0 1 1]);  
%  exportfig(gcf,['raw_data_' num2str(printver)]);
% exportfig(gcf,[foldername 'nlmod_mu']);
   yticks([0.4:0.1:0.6]);
    set(gca,'YLim',[0.35 0.65])
exportfig(gcf,[foldername 'nlmod_mu'],'color','cmyk','resolution',1600);

 
%% RATES
 
figure; % Endogenous policy rates
plot(x_vec ,(Rm_us_vec.^(freq)-1)*100,'b', 'linewidth',2.5); hold on;
plot(x_vec,(Rm_eu_vec.^(freq)-1)*100, 'linewidth',2.5,'LineStyle','--','color',col_eu);
%splot1(x_vec,(Rm_us_vec.^(freq)-1)*100*100); hold on;
%splot2(x_vec,(Rm_eu_vec.^(freq)-1)*100*100); 
xlim([x_vec(1) x_vec(end)]);
xlab=xlabel(x_lab,'Fontsize',FSize)
  % ylabel('%','interpreter','latex');
%formataxis(gca);
h=legend('$R^{m,us}$','$R^{m,eu}$','color','none','box','off','location','northeast','FontSize',20);
 
set(h,    'Position',[0.684931195350353 0.560912163246958 0.203210158575149 0.163095238095238])
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', FSize);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
ax = gca;
% Specify the number of desired ticks on the x-axis
% desiredNumXTicks = num_Tick; % Change this to the number of ticks you want
% xLimits = get(ax, 'XLim');
% customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
% set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
% xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
% set(gcf, 'PaperUnits', 'normalized');
% orient landscape;
% set(gcf, 'PaperPosition', [0 0 1 1]);  
%   set(gca,'Fontsize',FS_plot)
%  exportfig(gcf,['raw_data_' num2str(printver)]);
  set(gca,'Fontsize',FS_plot)
  xticks(xticks_def)
set(gca,'XLim',xlim_def)
 yticks([-0.30:0.15:0.15])
   set(gca,'YLim',[-0.30 0.15])
exportfig(gcf,[foldername 'nlmod_Rm'],'color','cmyk','resolution',1600);


figure;
plot(x_vec ,(Rd_us_vec.^(freq)-1)*100,'b', 'linewidth',2.5); hold on;
plot(x_vec,(Rd_eu_vec.^(freq)-1)*100, 'linewidth',2.5,'LineStyle','--','color',col_eu);
% 
% splot1(x_vec,(Rd_us_vec.^(freq)-1)*100*100); hold on;
% splot2(x_vec,(Rd_eu_vec.^(freq)-1)*100*100); 
xlim([x_vec(1) x_vec(end)]);
xlab=xlabel(x_lab,'Fontsize',FSize)
% ylabel('%','interpreter','latex');
%formataxis(gca);
h=legend('$R^{d,us}$','$R^{d,eu}$','color','none','box','off','location','southwest','FontSize',20);
set(h,    'Position',[0.684931195350353 0.560912163246958 0.203210158575149 0.163095238095238])
set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', FSize);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
ax = gca;
% Specify the number of desired ticks on the x-axis
% desiredNumXTicks = num_Tick; % Change this to the number of ticks you want
% xLimits = get(ax, 'XLim');
% customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
% set(ax, 'XTick', customXTicks);
% yLimits = get(ax, 'YLim');
% customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
% set(ax, 'YTick', customYTicks);
% desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
% xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
% set(gcf, 'PaperUnits', 'normalized');
% orient landscape;
%   set(gca,'Fontsize',FS_plot)
% set(gcf, 'PaperPosition', [0 0 1 1]);  
%  exportfig(gcf,['raw_data_' num2str(printver)]);
  set(gca,'Fontsize',FS_plot)
  xticks(xticks_def)
set(gca,'XLim',xlim_def)
  set(gca,'YLim',[-2.5 0.2])
exportfig(gcf,[foldername 'nlmod_Rd'],'color','cmyk','resolution',1600);

%% exchange rate


figure;
plot(x_vec,e_euus_vec,'b', 'linewidth',2.5); hold on;
xlab=xlabel(x_lab,'Fontsize',FSize)
%ylabel('bps','interpreter','latex');
%formataxis(gca);
xlim([x_vec(1) x_vec(end)]);
% h=legend('$R^{m,us}$','$R^{m,eu}$','color','none','box','off','location','northeast','FontSize',20);
set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', FSize);
% title('Exchange Rate','interpreter','latex','Fontsize',15);
ax = gca;
% Specify the number of desired ticks on the x-axis
% desiredNumXTicks = num_Tick; % Change this to the number of ticks you want
% xLimits = get(ax, 'XLim');
% customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
% set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
 % xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
  set(gca,'Fontsize',FS_plot)
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
% set(gcf, 'PaperUnits', 'normalized');
% orient landscape;
% set(gcf, 'PaperPosition', [0 0 1 1]);  
%  exportfig(gcf,['raw_data_' num2str(printver)]);
xticks(xticks_def)
set(gca,'XLim',xlim_def)
   yticks([2.1:0.1:2.3]);
    % ytickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
exportfig(gcf,[foldername 'nlmod_e'],'color','cmyk','resolution',1600);


 
