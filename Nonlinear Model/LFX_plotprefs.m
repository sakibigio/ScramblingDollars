%% [XII] Main Plots
% plot coordinates
path_g       = [cd '\Figs'];
imprime      = @(x) print( gcf, '-depsc2', [path_g filesep x]);
imprpdf      = @(x) eps2pdf( [path_g filesep x '.eps']);
printsb      = @(x) print(gcf,'-dpdf','-fillpage',[path_g filesep x]);
formataxis   = @(x) set(x, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 22, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
formatlegend = @(x) set(legend, 'Location', x, 'Orientation', 'Vertical', 'Box', 'Off', 'Fontsize', 22, 'Fontangle', 'normal');
label_x      = @(x) xlabel(x,'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 56,'Interpreter','latex');
label_y      = @(x) ylabel(x,'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 56,'interpreter','latex');
label_z      = @(x) zlabel(x,'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 56,'interpreter','latex');
otitle       = @(x) title(x, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 56,'interpreter','latex');


% Plot Preferences
color1   =  [0.4 0.2 0.8];
color2   =  [0.2 0.6 0.8];
color3   =  [0.1 0.4 0.4];
color4   =  [0.8 0.2 0.8];
color_mat=  [color1; color2; color3; color4];
color_base=[0.1 0.1 0.6]; 
color_base2=[0.6 0.1 0.1];
color_rss=  [0.6 0.1 0.1];
color_rss_val=  [0.6 0.6 0.1];
linewidthbase=6;
linewidthbase2=2;
markersize   = 120;

% Plot Lengths
Periods = (0:T-1);
pre_rss=-2	     ;
t_plotmax = 9   ;

% Plot Specifics
splot1=@(t,xt) plot(t, xt,'LineWidth',linewidthbase,'Color',color_base);
splot2=@(t,xt) plot(t, xt,'-.','LineWidth',linewidthbase,'Color',color_base2);
splot3=@(t,xt) plot(t, xt,'--','LineWidth',linewidthbase2,'Color',color3);
splot4=@(t,xt) plot(t, xt,':','LineWidth',linewidthbase2,'Color',color4);