% plot_paper_figures.m
% Extracted from main_LFX.m (lines 1788–2490, 2696–2834) on 2026-02-28
% Called by: main_LFX.m
%
% Contents:
%   - Policy function figures (mu, Rm, exchange rate, BP, CIP)
%   - Intensity plots
%   - Single regime plots (scrambling episodes)
%   - Phase diagram-like plots
%
% Figures saved to foldername (Overleaf quantfigs directory)

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



% Begin Figures
% figure('Name','Liquidity Ratios','NumberTitle','off'); % Liquidity Ratios
% splot1(x_vec,mu_us_vec); hold on;
% splot2(x_vec,mu_eu_vec);
% xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% % ylabel('$liquidity$','interpreter','latex');
% formataxis(gca);
% h=legend('$\mu^{us}$','$\mu^{eu}$','color','none','box','off','location','southwest','FontSize',20);
% set(h,'interpreter','latex');
% %scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
% %ylim(mean(mu_us_vec)+[-scale,scale]);
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% % set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% % title('Liquidity','interpreter','latex','Fontsize',15);
% ax = gca;
% % Specify the number of desired ticks on the x-axis
% desiredNumXTicks = 5; % Change this to the number of ticks you want
% xLimits = get(ax, 'XLim');
% customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
% set(ax, 'XTick', customXTicks);
% yLimits = get(ax, 'YLim');
% customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
% set(ax, 'YTick', customYTicks);
% desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% % Format the x-axis tick labels
% xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% % title('Real Policy Rates','interpreter','latex','Fontsize',15);
% set(gcf, 'PaperUnits', 'normalized');
% orient landscape;
% set(gcf, 'PaperPosition', [0 0 1 1]);  
% %  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
% print(gcf,'-dpdf',[foldername 'nlmod_mu']);

% Limit Graph

%limsigmas=[sigma_us_vec(1) sigma_us_vec(1)+0.6*(sigma_us_vec(end)-sigma_us_vec(1))];
limsigmas=[sigma_us_vec(1) sigma_us_vec(end-10)];


figure('Name','US Liquidity Ratio','NumberTitle','off'); % Liquidity Ratios
splot1(sigma_us_vec(index1),mu_us_vec(index1)); hold on;
splot2(sigma_us_vec(index2),mu_us_vec(index2));
xlim(limsigmas);
xlabel(x_lab,'interpreter','latex');
% ylabel('$liquidity$','interpreter','latex');
formataxis(gca);
% h=legend('$\mu^{us}$','$\mu^{eu}$','color','none','box','off','location','southwest','FontSize',20);
% set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Liquidity','interpreter','latex','Fontsize',15);
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 3; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  
%  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
if printit==1
    exportgraphics(gcf, fullfile(foldername, 'nlmod_mu_us.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

figure('Name','EU Liquidity Ratio','NumberTitle','off'); % Liquidity Ratios
splot1(sigma_us_vec(index1),mu_eu_vec(index1)); hold on;
splot2(sigma_us_vec(index2),mu_eu_vec(index2));
% xlim([sigma_us_vec(1) sigma_us_vec(end)]);
xlim(limsigmas);
xlabel(x_lab,'interpreter','latex');
% ylabel('$liquidity$','interpreter','latex');
formataxis(gca);
h=legend('$\mu^{us}$','$\mu^{eu}$','color','none','box','off','location','southwest','FontSize',20);
set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Liquidity','interpreter','latex','Fontsize',15);
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 5; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  
%  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
if printit==1
    exportgraphics(gcf, fullfile(foldername, 'nlmod_mu.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

% figure('Name','Real Rates Dollar','NumberTitle','off'); % Endogenous policy rates
% splot1(x_vec(index1),(Rm_us_vec(index1).^(freq)-1)*rate_scale); hold on;
% splot2(x_vec(index1),(Rm_us_vec(index1).^(freq)-1)*rate_scale); 
% xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('bps','interpreter','latex');
% formataxis(gca);
% h=legend('$R^{m,us}$','$R^{m,eu}$','color','none','box','off','location','northeast','FontSize',20);
% set(h,'interpreter','latex');
% %scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
% %ylim(mean(mu_us_vec)+[-scale,scale]);
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% % set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% % title('Real Policy Rates','interpreter','latex','Fontsize',15);
% ax = gca;
% % Specify the number of desired ticks on the x-axis
% desiredNumXTicks = 5; % Change this to the number of ticks you want
% xLimits = get(ax, 'XLim');
% customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
% set(ax, 'XTick', customXTicks);
% yLimits = get(ax, 'YLim');
% customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
% set(ax, 'YTick', customYTicks);
% desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% % Format the x-axis tick labels
% xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% % title('Real Policy Rates','interpreter','latex','Fontsize',15);
% set(gcf, 'PaperUnits', 'normalized');
% orient landscape;
% set(gcf, 'PaperPosition', [0 0 1 1]);  
% %  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
% print(gcf,'-dpdf',[foldername 'nlmod_Rm']);

figure('Name','Real Dollar Rate','NumberTitle','off'); % Endogenous policy rates
splot1(sigma_us_vec(index1),(Rm_us_vec(index1).^(freq)-1)*rate_scale); hold on;
splot2(sigma_us_vec(index2),(Rm_us_vec(index2).^(freq)-1)*rate_scale); 
% xlim([sigma_us_vec(1) sigma_us_vec(end)]);
xlim(limsigmas);
xlabel(x_lab,'interpreter','latex');
ylabel('bps','interpreter','latex');
formataxis(gca);
h=legend('Normal Regime','Scrambling Regime','color','none','box','off','location','northeast','FontSize',20);
set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 5; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  
%  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
if printit==1
    exportgraphics(gcf, fullfile(foldername, 'nlmod_Rm_us.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

figure('Name','Dollar Inflation','NumberTitle','off'); % Endogenous policy rates
splot1(sigma_us_vec(index1),(pi_us_vec(index1).^(freq)-1)*rate_scale); hold on;
splot2(sigma_us_vec(index2),(pi_us_vec(index2).^(freq)-1)*rate_scale); 
% xlim([sigma_us_vec(1) sigma_us_vec(end)]);
xlim(limsigmas);
xlabel(x_lab,'interpreter','latex');
ylabel('bps','interpreter','latex');
formataxis(gca);
h=legend('Normal Regime','Scrambling Regime','color','none','box','off','location','northeast','FontSize',20);
set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 5; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  
%  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
if printit==1
    exportgraphics(gcf, fullfile(foldername, 'nlmod_Rm_us.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

figure('Name','Real Euro Rate','NumberTitle','off'); % Endogenous policy rates
splot1(sigma_us_vec(index1),(Rm_eu_vec(index1).^(freq)-1)*rate_scale); hold on;
xlim([sigma_us_vec(1) sigma_us_vec(end)]);
splot2(sigma_us_vec(index2),(Rm_eu_vec(index2).^(freq)-1)*rate_scale); 
% xlim([sigma_us_vec(1) sigma_us_vec(end)]);
xlim(limsigmas);
xlabel(x_lab,'interpreter','latex');
ylabel('bps','interpreter','latex');
formataxis(gca);
h=legend('Normal Regime','Scrambling Regime','color','none','box','off','location','northeast','FontSize',20);
set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 5; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  
%  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
if printit==1
    exportgraphics(gcf, fullfile(foldername, 'nlmod_Rm_eu.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

figure('Name','all rates','NumberTitle','off'); % Endogenous policy rates
splot1(sigma_us_vec(index1),(Rm_us_vec(index1).^(freq)-1)*rate_scale); hold on;
splot2(sigma_us_vec(index2),(Rm_us_vec(index2).^(freq)-1)*rate_scale); 
splot3(sigma_us_vec(index1),(Rm_eu_vec(index1).^(freq)-1)*rate_scale); hold on;
splot4(sigma_us_vec(index2),(Rm_eu_vec(index2).^(freq)-1)*rate_scale); 
xlim([sigma_us_vec(1) sigma_us_vec(end)]);
xlabel(x_lab,'interpreter','latex');
ylabel('bps','interpreter','latex');
formataxis(gca);
h=legend('Normal Regime','Scrambling Regime','color','none','box','off','location','northeast','FontSize',20);
set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 5; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  


figure('Name','FX','NumberTitle','off');
e_mean=mean(e_euus_vec);
splot1(sigma_us_vec(index1),e_euus_vec(index1)/e_mean); hold on;
splot2(sigma_us_vec(index2),e_euus_vec(index2)/e_mean); 
xlabel(x_lab,'interpreter','latex');
formataxis(gca);
% xlim([sigma_us_vec(1) sigma_us_vec(end)]);
xlim(limsigmas);
h=legend('Normal Regime','Scrambling Regime','color','none','box','off','location','east','FontSize',20,'AutoUpdate', 'off');
set(h,'interpreter','latex');
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 3; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  
%  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
if printit==1
    exportgraphics(gcf, fullfile(foldername, 'nlmod_e.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end
area(sigma_us_vec(index1),yLimits(1)+(yLimits(2)-yLimits(1))*invp1/max([invp1; invp2]),'FaceAlpha',0.5,'FaceColor',color_base, 'DisplayName', 'Normal Regime'); hold on;
area(sigma_us_vec(index2),yLimits(1)+(yLimits(2)-yLimits(1))*invp2/max([invp1; invp2]),'FaceAlpha',0.5,'FaceColor',color_base2); hold on;
ylim([yLimits(1) yLimits(2)]);
xlim(limsigmas);
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  
if printit==1
    exportgraphics(gcf, fullfile(foldername, 'nlmod_FX.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

figure('Name','Invariants','NumberTitle','off');
area(sigma_us_vec(index1),invp1/max([invp1; invp2]),'FaceAlpha',0.5,'FaceColor',color_base, 'DisplayName', 'Normal Regime'); hold on;
area(sigma_us_vec(index2),invp2/max([invp1; invp2]),'FaceAlpha',0.5,'FaceColor',color_base2); hold on;
xlim(limsigmas);
xlabel(x_lab,'interpreter','latex');
formataxis(gca);
h=legend('Normal Regime','Scrambling Regime','color','none','box','off','location','northeast','FontSize',20,'AutoUpdate', 'off');
set(h,'interpreter','latex');
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 3; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
ylim([yLimits(1) yLimits(2)]);
if printit==1
    set(gcf, 'PaperUnits', 'normalized');
    orient landscape;
    set(gcf, 'PaperPosition', [0 0 1 1]);  
    exportgraphics(gcf, fullfile(foldername, 'nlmod_invariants.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

% figure('Name','Deposit Rates','NumberTitle','off');
% splot1(sigma_us_vec,(Rd_us_vec.^(freq)-1)*rate_scale); hold on;
% splot2(sigma_us_vec,(Rd_eu_vec.^(freq)-1)*rate_scale); 
% xlim([sigma_us_vec(1) sigma_us_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('bps','interpreter','latex');
% formataxis(gca);
% h=legend('$R^{d,us}$','$R^{d,eu}$','color','none','box','off','location','southwest','FontSize',20);
% set(h,'interpreter','latex');
% %scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
% %ylim(mean(mu_us_vec)+[-scale,scale]);
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% % set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% % title('Real Policy Rates','interpreter','latex','Fontsize',15);
% ax = gca;
% % Specify the number of desired ticks on the x-axis
% desiredNumXTicks = 5; % Change this to the number of ticks you want
% xLimits = get(ax, 'XLim');
% customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
% set(ax, 'XTick', customXTicks);
% yLimits = get(ax, 'YLim');
% customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
% set(ax, 'YTick', customYTicks);
% desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% % Format the x-axis tick labels
% xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% % title('Real Policy Rates','interpreter','latex','Fontsize',15);
% set(gcf, 'PaperUnits', 'normalized');
% orient landscape;
% set(gcf, 'PaperPosition', [0 0 1 1]);  
% %  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
% print(gcf,'-dpdf',[foldername 'nlmod_Rd']);

figure('Name','Bond Premium','NumberTitle','off')
splot1(sigma_us_vec(index1),(Rb_us_vec(index1).^(freq)-Rm_us_vec(index1).^(freq))*rate_scale); hold on;
xlim([sigma_us_vec(1) sigma_us_vec(end)]);
splot2(sigma_us_vec(index2),(Rb_us_vec(index2).^(freq)-Rm_us_vec(index2).^(freq))*rate_scale);
xlim([sigma_us_vec(1) sigma_us_vec(end)]);
xlabel(x_lab,'interpreter','latex');
ylabel('bps','interpreter','latex');
formataxis(gca);
%h=legend('$R^{d,us}$','$R^{d,eu}$','color','none','box','off','location','northeast','FontSize',20);
% set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 3; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  
%  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
if printit==1
    exportgraphics(gcf, fullfile(foldername, 'nlmod_BP.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

figure('Name','Real Rate Differential','NumberTitle','off')
plot(sigma_us_vec(index1),0*Rm_eu_vec(index1), 'LineStyle', '--','Color','k','LineWidth',1); hold on;
aux1=(Rm_eu_vec(index1).^(freq)-Rm_us_vec(index1).^(freq))*rate_scale;
aux2=(Rm_eu_vec(index2).^(freq)-Rm_us_vec(index2).^(freq))*rate_scale;
splot1(sigma_us_vec(index1),aux1); hold on;
splot2(sigma_us_vec(index2),aux2);
% xlim([sigma_us_vec(1) sigma_us_vec(end)]);
xlim(limsigmas);
xlabel(x_lab,'interpreter','latex');
ylabel('bps','interpreter','latex');
formataxis(gca);
%h=legend('$R^{d,us}$','$R^{d,eu}$','color','none','box','off','location','northeast','FontSize',20);
% set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
% Get the current axes handle
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 3; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
if printit==1
    %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
    % title('Real Policy Rates','interpreter','latex','Fontsize',15);
    set(gcf, 'PaperUnits', 'normalized');
    orient landscape;
    set(gcf, 'PaperPosition', [0 0 1 1]);  
    %  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
    exportgraphics(gcf, fullfile(foldername, 'nlmod_UIP.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

figure('Name','Real Rate Differential shade','NumberTitle','off')
plot(sigma_us_vec(index1),0*Rm_eu_vec(index1), 'LineStyle', '--','Color','k','LineWidth',1); hold on;
aux1=(Rm_eu_vec(index1).^(freq)-Rm_us_vec(index1).^(freq))*rate_scale;
aux2=(Rm_eu_vec(index2).^(freq)-Rm_us_vec(index2).^(freq))*rate_scale;
diff=aux2-aux1;
for ii=1:length(index1)-1
    alphaVal1 = invp1(ii)/max(invp1);
    alphaVal2 = invp2(ii)/max(invp2);
    plot(sigma_us_vec(index1(ii:ii+1)),aux1(ii:ii+1),'Color',[0 0 1 alphaVal1],'LineWidth',2); hold on;
    plot(sigma_us_vec(index2(ii:ii+1)),aux2(ii:ii+1),'Color',[1 0 1 alphaVal2],'LineWidth',2);
end
%splot3(sigma_us_vec(index2),diff);
xlim([sigma_us_vec(1) sigma_us_vec(end)]);
xlabel(x_lab,'interpreter','latex');
ylabel('bps','interpreter','latex');
formataxis(gca);
%h=legend('$R^{d,us}$','$R^{d,eu}$','color','none','box','off','location','northeast','FontSize',20);
% set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
% Get the current axes handle
ax = gca;

%% Intensity Plots
figure('Name','Real Rate Differential shade','NumberTitle','off')
plot(sigma_us_vec(index1),0*e_euus_vec(index1)+1, 'LineStyle', '--','Color','k','LineWidth',1); hold on;
aux1=e_euus_vec(index1)/mean(e_euus_vec);
aux2=e_euus_vec(index2)/mean(e_euus_vec);
for ii=1:length(index1)-1
    alphaVal1 = sqrt(invp1(ii)/max(invp1));
    alphaVal2 = sqrt(invp2(ii)/max(invp2));
    plot(sigma_us_vec(index1(ii:ii+1)),aux1(ii:ii+1),'Color',[0 0 1 alphaVal1],'LineWidth',2); hold on;
    plot(sigma_us_vec(index2(ii:ii+1)),aux2(ii:ii+1),'Color',[1 0 1 alphaVal2],'LineWidth',2);
end
%splot3(sigma_us_vec(index2),diff);
xlim([sigma_us_vec(1) sigma_us_vec(end)]);
xlabel(x_lab,'interpreter','latex');
ylabel('bps','interpreter','latex');
formataxis(gca);
%h=legend('$R^{d,us}$','$R^{d,eu}$','color','none','box','off','location','northeast','FontSize',20);
% set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
% Get the current axes handle
ax = gca;

figure('Name','Real Rate Differential shade','NumberTitle','off')
plot(sigma_us_vec(index1),0*Rm_eu_vec(index1), 'LineStyle', '--','Color','k','LineWidth',1); hold on;
aux1=(Rm_eu_vec(index1).^(freq)-Rm_us_vec(index1).^(freq))*rate_scale;
aux2=(Rm_eu_vec(index2).^(freq)-Rm_us_vec(index2).^(freq))*rate_scale;
diff=aux2-aux1;
for ii=1:length(index1)-1
    alphaVal1 = sqrt(invp1(ii)/max(invp1));
    alphaVal2 = sqrt(invp2(ii)/max(invp2));
    plot(sigma_us_vec(index1(ii:ii+1)),aux1(ii:ii+1),'Color',[0 0 1 alphaVal1],'LineWidth',2); hold on;
    plot(sigma_us_vec(index2(ii:ii+1)),aux2(ii:ii+1),'Color',[1 0 1 alphaVal2],'LineWidth',2);
end
%splot3(sigma_us_vec(index2),diff);
xlim([sigma_us_vec(1) sigma_us_vec(end)]);
xlabel(x_lab,'interpreter','latex');
ylabel('bps','interpreter','latex');
formataxis(gca);
%h=legend('$R^{d,us}$','$R^{d,eu}$','color','none','box','off','location','northeast','FontSize',20);
% set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
% Get the current axes handle
ax = gca;


%%  Single Regime Plots
figure('Name','FX','NumberTitle','off');
% e_mean=mean(e_euus_vec);
% splot1(sigma_us_vec(index1),e_euus_vec(index1)/e_mean); hold on;
splot1(sigma_us_vec(index2),e_euus_vec(index2)/e_mean); 
xlabel(x_lab,'interpreter','latex');
formataxis(gca);
xlim([sigma_us_vec(1) sigma_us_vec(end)]);
%h=legend('Normal Regime','Scrambling Regime','color','none','box','off','location','east','FontSize',20,'AutoUpdate', 'off');
% set(h,'interpreter','latex');
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 3; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  
%  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
if printit==1
    exportgraphics(gcf, fullfile(foldername, 'nlmod_e_scramble.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end
area(sigma_us_vec(index1),yLimits(1)+(yLimits(2)-yLimits(1))*invp1/max([invp1; invp2]),'FaceAlpha',0.5,'FaceColor',color_base, 'DisplayName', 'Normal Regime'); hold on;
area(sigma_us_vec(index2),yLimits(1)+(yLimits(2)-yLimits(1))*invp2/max([invp1; invp2]),'FaceAlpha',0.5,'FaceColor',color_base2); hold on;
ylim([yLimits(1) yLimits(2)]);
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  
if printit==1
    exportgraphics(gcf, fullfile(foldername, 'nlmod_FX_scramble.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

figure('Name','Real Rate Differential','NumberTitle','off')
plot(sigma_us_vec(index1),0*Rm_eu_vec(index1), 'LineStyle', '--','Color','k','LineWidth',1); hold on;
%aux1=(Rm_eu_vec(index1).^(freq)-Rm_us_vec(index1).^(freq))*rate_scale;
aux2=(Rm_eu_vec(index2).^(freq)-Rm_us_vec(index2).^(freq))*rate_scale;
%diff=aux2-aux1;
splot1(sigma_us_vec(index1),aux2); hold on;
% splot2(sigma_us_vec(index2),aux2);
%splot3(sigma_us_vec(index2),diff);
xlim([sigma_us_vec(1) sigma_us_vec(end)]);
xlabel(x_lab,'interpreter','latex');
ylabel('bps','interpreter','latex');
formataxis(gca);
%h=legend('$R^{d,us}$','$R^{d,eu}$','color','none','box','off','location','northeast','FontSize',20);
% set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
% Get the current axes handle
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 3; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
if printit==1
    %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
    % title('Real Policy Rates','interpreter','latex','Fontsize',15);
    set(gcf, 'PaperUnits', 'normalized');
    orient landscape;
    set(gcf, 'PaperPosition', [0 0 1 1]);  
    %  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
    exportgraphics(gcf, fullfile(foldername, 'nlmod_UIP_scramble.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

figure('Name','US Liquidity Ratio','NumberTitle','off'); % Liquidity Ratios
% splot1(sigma_us_vec(index1),mu_us_vec(index1)); hold on;
splot1(sigma_us_vec(index2),mu_us_vec(index2));
xlim([sigma_us_vec(1) sigma_us_vec(end)]);
xlabel(x_lab,'interpreter','latex');
% ylabel('$liquidity$','interpreter','latex');
formataxis(gca);
% h=legend('$\mu^{us}$','$\mu^{eu}$','color','none','box','off','location','southwest','FontSize',20);
% set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Liquidity','interpreter','latex','Fontsize',15);
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 3; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  
%  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
if printit==1
    exportgraphics(gcf, fullfile(foldername, 'nlmod_mu_us_scramble.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end


%% Phase Diagram-Like Plots
% Interpolation
N_r1_1=1;
N_r2_1=60;
N_r1_2=20;
N_sim =N_r1_1+N_r2_1+N_r1_2;
lsigma_0=mu_sigma_us_r2+2;
lsigma_typath=zeros(N_sim,1);
lsigma_typath(1)=lsigma_0;
for ii=1:N_sim
    if ii<=N_r1_1
        lsigma_typath(ii+1)=mu_sigma_us_r1*(1-rho_sigma_us_r1)+rho_sigma_us_r1*lsigma_typath(ii,1);
    elseif ii>(N_r1_1+N_r2_1)
        lsigma_typath(ii+1)=mu_sigma_us_r1*(1-rho_sigma_us_r1)+rho_sigma_us_r1*lsigma_typath(ii,1);
    else
        lsigma_typath(ii+1)=mu_sigma_us_r2*(1-rho_sigma_us_r2)+rho_sigma_us_r2*lsigma_typath(ii,1);
    end
end
sigma_typath=exp(lsigma_typath);
e_euus_typath=0*sigma_typath;
for ii=1:N_sim
    if ii<=N_r1_1
        e_euus_typath(ii)=interp1(sigma_us_vec(index1),e_euus_vec(index1)/e_mean,sigma_typath(ii),'spline');
    elseif ii>N_r1_1+N_r2_1
        e_euus_typath(ii)=interp1(sigma_us_vec(index1),e_euus_vec(index1)/e_mean,sigma_typath(ii),'spline');
    else
        e_euus_typath(ii)=interp1(sigma_us_vec(index2),e_euus_vec(index2)/e_mean,sigma_typath(ii),'spline');
    end
end

figure('Name','FX Path','NumberTitle','off');
e_mean=mean(e_euus_vec);
splot1(sigma_us_vec(index1),e_euus_vec(index1)/e_mean); hold on;
splot2(sigma_us_vec(index2),e_euus_vec(index2)/e_mean); 
for ii=1:N_sim-1
   plot(sigma_typath(ii:ii+1),e_euus_typath(ii:ii+1),'Color',[0.5 0 0.5],'LineWidth',1,'Marker','<','MarkerFaceColor', [0.5 0 0.5], 'MarkerEdgeColor', [0.5 0 0.5], 'MarkerSize', 10); hold on;
   % annotation('arrow',[sigma_typath(ii:ii+1)],[e_euus_typath(ii:ii+1)]); hold on;
end
xlabel(x_lab,'interpreter','latex');
formataxis(gca);
xlim([sigma_us_vec(1) sigma_us_vec(end)]);
h=legend('Normal Regime','Scrambling Regime','color','none','box','off','location','east','FontSize',20,'AutoUpdate', 'off');
set(h,'interpreter','latex');
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 3; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
hold on;

% figure('Name','FX Typical Path','NumberTitle','off')
% lsigma_0=mu_sigma_us_r1;
% lsigma_typath=zeros(36,1);
% lsigma_typath(1)=lsigma_0;
% %mu_sigma_us_r1 ;
% %rho_sigma_us_r1;
% for ii=1:5
%     aux=mu_sigma_us_r1*(1-rho_sigma_us_r1)+rho_sigma_us_r1*lsigma_typath(ii,1);
%     [~,index]=min(abs(sigma_us_vec(1:N_sigma_us/2)-exp(aux)));
%     lsigma_typath(ii+1)=log(sigma_us_vec(index));
% end
% for ii=6:22
%     aux=mu_sigma_us_r2*(1-rho_sigma_us_r2)+rho_sigma_us_r2*lsigma_typath(ii,1);
%     [~,index]=min(abs(sigma_us_vec(N_sigma_us/2+1:end)-exp(aux)));
%     lsigma_typath(ii+1)=log(sigma_us_vec(N_sigma_us/2+index));
% end
% for ii=23:30
%     aux=mu_sigma_us_r1*(1-rho_sigma_us_r1)+rho_sigma_us_r1*lsigma_typath(ii,1);
%     [~,index]=min(abs(sigma_us_vec(1:N_sigma_us/2)-exp(aux)));
%     lsigma_typath(ii+1)=log(sigma_us_vec(index));
% end
% %splot3(sigma_us_vec(index2),diff);
% xlim([sigma_us_vec(1) sigma_us_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('bps','interpreter','latex');
% formataxis(gca);
% 
% figure('Name','FX Typical Path','NumberTitle','off')
% sigma_0=mu_sigma_us_r1;
% sigma_typath=zeros(30,1);
% index=zeros(30,1);
% sigma_typath(1)=sigma_0;
% aux=mu_sigma_us_r1*(1-rho_sigma_us_r1)+rho_sigma_us_r1*sigma_typath(1);
% [~,index0]=min(abs(sigma_us_vec(1:N_sigma_us/2)-exp(aux)));
% index(1)=index0;
% sigma_typath(1)=sigma_us_vec(index0);
% %mu_sigma_us_r1 ;
% %rho_sigma_us_r1;
% for ii=1:5
%     [~,index_aux]=max(Zprob_sigma_us(index(ii),1:N_sigma_us/2));
%     index(ii+1)=index_aux;
%     sigma_typath(ii+1)=sigma_us_vec(index(ii+1));
% end
% [~,index0]=min(abs(sigma_us_vec(N_sigma_us/2+1:end)-(sigma_typath(ii+1))));
% index(ii+1)=N_sigma_us/2+index0;
% sigma_typath(ii+1)=sigma_us_vec(index(ii+1));
% for ii=6:22
%     [~,index_aux]=max(Zprob_sigma_us(index(ii),N_sigma_us/2+1:end));
%     index(ii+1)=N_sigma_us/2+index_aux;
%     sigma_typath(ii+1)=sigma_us_vec(index(ii+1));
% end
% [~,index0]=min(abs(sigma_us_vec(1:N_sigma_us/2)-(sigma_typath(ii+1))));
% index(ii+1)=index0;
% sigma_typath(ii+1)=sigma_us_vec(index(ii+1));
% for ii=23:30
%     [~,index_aux]=max(Zprob_sigma_us(index(ii),1:N_sigma_us/2));
%     index(ii+1)=index_aux;
%     sigma_typath(ii+1)=sigma_us_vec(index(ii+1));
% end
% figure
% for ii=1:5
%     line(sigma_us_vec(index1(ii:ii+1)),e_euus_vec(index1(ii:ii+1)),'Color',[1 0 0],'LineWidth',2,'Marker','<'); hold on;
% end
% for ii=6
%     line(sigma_us_vec(index1(ii)),e_euus_vec(index2(ii:ii+1)),'Color',[0 0 1],'LineWidth',2,'Marker','>'); hold on;
% end
% for ii=7:22
%     line(sigma_us_vec(index2(ii:ii+1)),e_euus_vec(index2(ii:ii+1)),'Color',[0 0 1],'LineWidth',2,'Marker','>'); hold on;
% end
% for ii=23:30
%     line(sigma_us_vec(index1(ii:ii+1)),e_euus_vec2(index(ii:ii+1)),'Color',[1 0 0],'LineWidth',2,'Marker','<'); hold on;
% end
% %splot3(sigma_us_vec(index2),diff);
% xlim([sigma_us_vec(1) sigma_us_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('bps','interpreter','latex');
% formataxis(gca);