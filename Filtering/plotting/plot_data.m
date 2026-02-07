%% Plot data
% Set output folder based on machine
[~, username] = system('whoami');
username = strtrim(username);
if strcmp(username, 'sakibigio')
    foldername = '/Users/sakibigio/Dropbox/Apps/Overleaf/ScramblingDollars_Revision_Restud/FigsTabs/';
elseif strcmp(username, 'sakiclaudia')
    foldername = '/Users/sakiclaudia/Library/CloudStorage/Dropbox/Apps/Overleaf/ScramblingDollarsLiquidity_NewVersion/FigsTabs/';
else
    error('Unknown user: %s. Please set foldername manually.', username);
end

load LFX_data3;
load exchange_rate_data;

%% Plot estimated sequences
formataxis   = @(x) set(x, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 28, 'Box','On','PlotBoxAspectRatio', [1 0.75 1]);
label_x      = @(x) xlabel(x,'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 28,'Interpreter','latex');
label_y      = @(x) ylabel(x,'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 28,'interpreter','latex');

%dates=datenum(2001,1:192,1);
%datesperiod = 36:160;
% datesperiod = 1:192;

datesperiod = 1:234;
dates=datenum(2001,1:234,1);
abs_scale=12e4;

curlist = {'au','ca','jp','nz','no','sw','ch','uk'};
conlist = {'AUD','CAD','JPY','NZD','NOK','SWK','CHF','GBP'};
CURRlist = {'EUR','AUD','CAD','JPY','NZD','NOK','SWK','CHF','GBP'};

%% Separate subfigures
figure('Name','Dollar Euro FX','NumberTitle','off')   
plot(dates(datesperiod),exp(-inv_e(datesperiod)),'LineWidth',2);
datetick('x','yyyy-mm','keeplimits');
label_x('Time (Year-Month)');
formataxis(gca);
xlabel('year','interpreter','latex');
ylabel('level','interpreter','latex');
grid on; axis tight;
% title('$e $ (EUR/USD)','interpreter','latex','Fontsize',20);
if printit==1
    % orient landscape;
    set(gcf, 'PaperUnits', 'normalized');
    orient landscape;
    set(gcf, 'PaperPosition', [0 0 1 1]);  
    %  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
    % orient landscape;
    print(gcf,'-dpdf',[foldername 'raw_data_FX']);
end

figure('Name','Liquidity Ratios','NumberTitle','off')   
plot(dates(datesperiod),exp(mu_us(datesperiod)),'LineWidth',2); hold on;
plot(dates(datesperiod),exp(mu_eu(datesperiod)),'LineWidth',2); hold off;
h=legend('US','EU','color','none','box','off','location','northwest');
set(h,'interpreter','latex');
datetick('x','yyyy-mm','keeplimits');
formataxis(gca);
xlabel('year','interpreter','latex');
ylabel('level','interpreter','latex');
grid on; axis tight;
xlabel('year','interpreter','latex');
ylabel('level','interpreter','latex');
% title('$\mu^{c}$','interpreter','latex','Fontsize',20);
if printit==1
    set(gcf, 'PaperUnits', 'normalized');
    orient landscape;
    set(gcf, 'PaperPosition', [0 0 1 1]);  
    %  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
    % orient landscape;
    print(gcf,'-dpdf',[foldername 'raw_data_mu']);
end

figure('Name','Rate Differentials','NumberTitle','off')   
plot(dates(datesperiod),(exp(im_eu(datesperiod))-1)*abs_scale,'LineWidth',2); hold on;
plot(dates(datesperiod),(exp(im_us(datesperiod))-1)*abs_scale,'LineWidth',2);
datetick('x','yyyy-mm','keeplimits'); formataxis(gca);
xlabel('year','interpreter','latex');
ylabel('annual bps','interpreter','latex');
grid on; axis tight;
% title('$i_m^{c}$','interpreter','latex','Fontsize',20);
h=legend('US','EU','color','none','box','off','location','northeast');
set(h,'interpreter','latex');
if printit==1
    set(gcf, 'PaperUnits', 'normalized');
    orient landscape;
    set(gcf, 'PaperPosition', [0 0 1 1]);  
    %  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
    print(gcf,'-dpdf',[foldername 'raw_data_im']);
end

figure('Name','Money Supplies','NumberTitle','off')   
plot(dates(datesperiod),exp(M_us(datesperiod)),'LineWidth',2); hold on;
plot(dates(datesperiod),exp(M_eu(datesperiod)),'LineWidth',2); hold off;
xlabel('year','interpreter','latex');
ylabel('level','interpreter','latex');
h=legend('US','EU','color','none','box','off','location','northwest');
set(h,'interpreter','latex');
datetick('x','yyyy-mm','keeplimits'); formataxis(gca);
grid on; axis tight;
% title('$M^{c}$','interpreter','latex','Fontsize',20);
if printit==1
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  
%  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
print(gcf,'-dpdf',[foldername 'raw_data_M']);
end

figure('Name','Bond and CIP premia (US)','NumberTitle','off')   
% subplot(3,2,5);
plot(dates(datesperiod),Rb_Rm(datesperiod)*abs_scale,'LineWidth',2); hold on;
plot(dates(datesperiod),cip(datesperiod)*abs_scale,'LineWidth',2); hold off;
h=legend('$\mathcal{BP}$','$\mathcal{CIP}$','color','none','box','off','location','northwest');
set(h,'interpreter','latex');
datetick('x','yyyy-mm','keeplimits'); formataxis(gca);
grid on; axis tight;
xlabel('year','interpreter','latex');
ylabel('annual bps','interpreter','latex');
%title('$R_b^*-R_m^*$','interpreter','latex','Fontsize',15);
%title('CIP','interpreter','latex','Fontsize',15);
% title('$\mathcal{BP},\mathcal{CIP}$','interpreter','latex','Fontsize',20);
if printit==1
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  
%  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
% 
print(gcf,'-dpdf',[foldername 'raw_data_spreads']);
end

figure('Name','Interbank Spreads (US)','NumberTitle','off')   
plot(dates(datesperiod),Chi_D_US(datesperiod)*abs_scale,'LineWidth',2);
datetick('x','yyyy-mm','keeplimits'); formataxis(gca);
grid on; axis tight;
xlabel('year','interpreter','latex');
ylabel('annual bps','interpreter','latex');
% title('Interbank Spread','interpreter','latex','Fontsize',20);
% h=legend('$\mathcal{OIS}$','$R^f Spread $','color','none','box','off','location','northwest');
%set(h,'interpreter','latex');
%h=suptitle('Data');
%set(h,'interpreter','latex','Fontsize',20);
if printit==1
    set(gcf, 'PaperUnits', 'normalized');
    orient landscape;
    set(gcf, 'PaperPosition', [0 0 1 1]);  
    %  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
    % orient landscape;
    print(gcf,'-dpdf',[foldername 'raw_data_CIP']);
    % if printit==1
    %     % orient landscape;
    %     set(gcf, 'PaperUnits', 'normalized');
    %     set(gcf, 'PaperPosition', [0 0 1 1]);  
    %   %  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
    %     print(gcf,'-dpdf',[foldername 'raw_data']);
end

% % Plot data variables separately
% figure
% plot(dates(datesperiod),exp(mu_us(datesperiod)),'LineWidth',2);
% datetick('x','yyyy-mm','keeplimits');
% grid on; axis tight;
% label_x('Time (Year-Month)');
% formataxis(gca);
% if printit==1
%     orient landscape;
%     set(gcf, 'PaperUnits', 'normalized');
%     set(gcf, 'PaperPosition', [0 0 1 1]);  
%     %print(gcf,'-dpdf',['fig_data_mbs_mu_us' num2str(printver)]);
%     h = title('Data: $\mu_{us}$');
%     set(h,'interpreter','latex','fontsize',20);
% end
% 
% figure
% plot(dates(datesperiod),exp(-inv_e(datesperiod)),'LineWidth',2);
% datetick('x','yyyy-mm','keeplimits');
% grid on; axis tight;
% label_x('Time (Year-Month)');
% formataxis(gca);
% if printit==1
%     orient landscape;
%     set(gcf, 'PaperUnits', 'normalized');
%     set(gcf, 'PaperPosition', [0 0 1 1]);  
%     %print(gcf,'-dpdf',['fig_data_mbs_euus' num2str(printver)]);
%     h = title('Data: log(EUR/USD)');
%     set(h,'interpreter','latex','fontsize',20);
% end
% 
% figure
% plot(dates(datesperiod),Rb_Rm(datesperiod)*12*1e4,'LineWidth',2);
% datetick('x','yyyy-mm','keeplimits');
% grid on; axis tight;
% label_x('Time (Year-Month)');
% label_y('Bps');
% formataxis(gca);
% if printit==1
%     orient landscape;
%     set(gcf, 'PaperUnits', 'normalized');
%     set(gcf, 'PaperPosition', [0 0 1 1]);  
%     %print(gcf,'-dpdf',['fig_data_mbs_Rb_Rm' num2str(printver)]);
%     h = title('Data: $R_b^{us}-R_m^{us}$');
%     set(h,'interpreter','latex','fontsize',20);
% end
% 
% figure
% plot(dates(datesperiod),cip(datesperiod)*12*1e4,'LineWidth',2);
% datetick('x','yyyy-mm','keeplimits');
% grid on; axis tight;
% label_x('Time (Year-Month)');
% label_y('Bps');
% formataxis(gca);
% if printit==1
%     orient landscape;
%     set(gcf, 'PaperUnits', 'normalized');
%     set(gcf, 'PaperPosition', [0 0 1 1]);  
%     %print(gcf,'-dpdf',['fig_data_mbs_CIP' num2str(printver)]);
%     h = title('Data: CIP');
%     set(h,'interpreter','latex','fontsize',20);
% end
% 
% figure
% plot(dates(datesperiod(12:end)),ois(datesperiod(12:end)-11)*12*1e4,'LineWidth',2);
% datetick('x','yyyy-mm','keeplimits');
% grid on; axis tight;
% label_x('Time (Year-Month)');
% label_y('Bps');
% formataxis(gca);
% if printit==1
%     orient landscape;
%     set(gcf, 'PaperUnits', 'normalized');
%     set(gcf, 'PaperPosition', [0 0 1 1]);  
%     %print(gcf,'-dpdf',['fig_data_mbs_OIS' num2str(printver)]);
%     h = title('Data: OIS');
%     set(h,'interpreter','latex','fontsize',20);
% end
% 
% figure
% plot(dates(datesperiod),exp(mu_eu(datesperiod)),'LineWidth',2);
% datetick('x','yyyy-mm','keeplimits');
% grid on; axis tight;
% label_x('Time (Year-Month)');
% formataxis(gca);
% if printit==1
%     orient landscape;
%     set(gcf, 'PaperUnits', 'normalized');
%     set(gcf, 'PaperPosition', [0 0 1 1]);  
%     %print(gcf,'-dpdf',['fig_data_mbs_mu_eu' num2str(printver)]);
%     h = title('Data: $\mu_{eu}$');
%     set(h,'interpreter','latex','fontsize',20);
% end
% figure
% plot(dates(datesperiod),(Chi_D_US(datesperiod))*abs_scale,'LineWidth',2);
% datetick('x','yy','keeplimits');
% grid on; axis tight;
% label_x('Time (Year-Month)');
% formataxis(gca);
% if printit==1
%     orient landscape;
%     set(gcf, 'PaperUnits', 'normalized');
%     set(gcf, 'PaperPosition', [0 0 1 1]);  
%     %print(gcf,'-dpdf',['fig_data_mbs_Chi_d_us' num2str(printver)]);
%     h = title('$\chi^-_{us}-\chi^+_{us}$','interpreter','latex','Fontsize',15);
%     set(h,'interpreter','latex','fontsize',20);
% end

% Plot estimated M_us, sigma_us, mu_us, and inv_e together


%% Interest Rate differentials
% Simulated exchange rates vs data (in levels)
figure
abp_scale=12*1e4;
subplot(3,3,1);plot(dates(datesperiod),(im_eu(datesperiod)-im_us(datesperiod))*abp_scale,'LineWidth',2);hold on;
datetick('x','yy');
grid on; axis tight;
label_x('Time (Year)');
formataxis(gca);
h=legend('','Data');
set(h,'interpreter','latex','location','Northeast','Fontsize',10);
title('EUR/USD','interpreter','latex','Fontsize',15);
for j=1:length(curlist)
    aux = eval(['im_' curlist{j} '(datesperiod)'])-im_us((datesperiod));
    %temp = eval(['ln_' curlist{j} '_us_t((datesperiod(1)))'])-eval(['-oo_.UpdatedVariables.inv_e_' curlist{j} '(datesperiod(1))']);
    subplot(3,3,j+1);
    plot(dates(datesperiod),aux*abp_scale,'LineWidth',2);hold on;
    %subplot(3,3,j+1);plot(eval(['-inv_e_' curlist{j} '(end-200:end)']),'LineWidth',2);
    datetick('x','yy');
    grid on; axis tight;
    label_x('Time (Year)');
    formataxis(gca);
    title([conlist{j} '/USD'],'interpreter','latex','Fontsize',15);
end
if printit==1
    orient landscape;
    set(gcf, 'PaperUnits', 'normalized');
    set(gcf, 'PaperPosition', [0 0 1 1]);
    print(gcf,'-dpdf',['fig_exchange_rate_path_mbs' num2str(printver)]);
end

% Ted Spreads
figure
abp_scale=12*1e4;
subplot(3,3,1); plot(dates(datesperiod),TED_s_eu_t(datesperiod)*abp_scale,'LineWidth',2);hold on;
plot(dates(datesperiod),TED_s_us_t(datesperiod)*abp_scale,'LineWidth',2);hold on;
datetick('x','yy');
grid on; axis tight;
label_x('Time (Year)');
formataxis(gca);
h=legend('Country','US');
set(h,'interpreter','latex','location','Northeast','Fontsize',10);
title('TED differentials','interpreter','latex','Fontsize',15);
for j=1:length(curlist)
    aux = eval(['TED_s_' curlist{j} '_t(datesperiod)']);
    %temp = eval(['ln_' curlist{j} '_us_t((datesperiod(1)))'])-eval(['-oo_.UpdatedVariables.inv_e_' curlist{j} '(datesperiod(1))']);
    subplot(3,3,j+1);
    plot(dates(datesperiod),aux*abp_scale,'LineWidth',2);hold on;
    plot(dates(datesperiod),TED_s_us_t(datesperiod)*abp_scale,'LineWidth',2);
    %subplot(3,3,j+1);plot(eval(['-inv_e_' curlist{j} '(end-200:end)']),'LineWidth',2);
    datetick('x','yy');
    grid on; axis tight;
    label_x('Time (Year)');
    formataxis(gca);
    title([conlist{j} ' Ted Spread'],'interpreter','latex','Fontsize',15);
end
if printit==1
    orient landscape;
    set(gcf, 'PaperUnits', 'normalized');
    set(gcf, 'PaperPosition', [0 0 1 1]);
    print(gcf,'-dpdf',['fig_exchange_rate_path_mbs' num2str(printver)]);
end

%% Non-Negativity test
% We are using that:
% im-im*+chi_us=chi_eu+rp
figure
abp_scale=12*1e4;
aux=(-(im_eu(datesperiod)-im_us(datesperiod))+Rb_Rm(datesperiod))*abp_scale;
rsp_aux=-(mean(aux)-mean(cip(datesperiod))*abp_scale);
m_aux=mean(aux)+rsp_aux;
subplot(3,3,1);
plot(dates(datesperiod),aux,'LineWidth',2); hold on;
plot(dates(datesperiod),m_aux+0*dates(datesperiod),'LineWidth',1,'Color','r'); hold off;
%subplot(3,3,1);plot(dates,-inv_e(end-191:end),'LineWidth',2);hold on;
datetick('x','yy');
grid on; axis tight;
label_x('Time (Year)');
formataxis(gca);
% h=legend('Model','Data');
% set(h,'interpreter','latex','location','Northeast','Fontsize',10);
title(['EUR/USD (' num2str(m_aux) ')'],'interpreter','latex','Fontsize',15);
for j=1:length(curlist)
    aux = (-(eval(['im_' curlist{j} '(datesperiod)'])-im_us(datesperiod))+Rb_Rm(datesperiod))*abp_scale;
    m_aux=mean(aux)+rsp_aux;
    %temp = eval(['ln_' curlist{j} '_us_t((datesperiod(1)))'])-eval(['-oo_.UpdatedVariables.inv_e_' curlist{j} '(datesperiod(1))']);
    subplot(3,3,j+1);
    plot(dates(datesperiod),aux,'LineWidth',2);hold on;
    plot(dates(datesperiod),m_aux+0*dates(datesperiod),'LineWidth',1,'Color','r'); hold off;
    %subplot(3,3,j+1);plot(eval(['-inv_e_' curlist{j} '(end-200:end)']),'LineWidth',2);
    datetick('x','yy');
    grid on; axis tight;
    label_x('Time (Year)');
    formataxis(gca);
    title([conlist{j} '/USD (' num2str(m_aux) ')'],'interpreter','latex','Fontsize',15);
end


if printit==1
    orient landscape;
    set(gcf, 'PaperUnits', 'normalized');
    set(gcf, 'PaperPosition', [0 0 1 1]);
    print(gcf,'-dpdf',['fig_exchange_rate_path_mbs' num2str(printver)]);
end

%% Single Figure Plots
% In a single figure
figure
subplot(3,3,1);plot(dates(datesperiod),exp(mu_us(datesperiod)),'LineWidth',2);
datetick('x','yy','keeplimits');
grid on; axis tight;
title('$\mu^*$','interpreter','latex','Fontsize',15);
subplot(3,3,2);plot(dates(datesperiod),exp(-inv_e(datesperiod)),'LineWidth',2);
datetick('x','yy','keeplimits');
grid on; axis tight;
title('$EUR/USD$','interpreter','latex','Fontsize',15);
subplot(3,3,3);plot(dates(datesperiod),Rb_Rm(datesperiod)*abs_scale,'LineWidth',2);
datetick('x','yy','keeplimits');
grid on; axis tight;
title('$R_b^*-R_m^*$','interpreter','latex','Fontsize',15);
subplot(3,3,4);plot(dates(datesperiod),exp(mu_eu(datesperiod)),'LineWidth',2);
datetick('x','yy','keeplimits');
grid on; axis tight;
title('$\mu$','interpreter','latex','Fontsize',15);
subplot(3,3,5);plot(dates(datesperiod),exp(im_eu(datesperiod)),'LineWidth',2);
datetick('x','yy','keeplimits');
grid on; axis tight;
title('$i_m$','interpreter','latex','Fontsize',15);
subplot(3,3,6);plot(dates(datesperiod),exp(im_us(datesperiod)),'LineWidth',2);
datetick('x','yy','keeplimits');
grid on; axis tight;
title('$i_m^*$','interpreter','latex','Fontsize',15);
subplot(3,3,7);plot(dates(datesperiod),(cip(datesperiod))*abs_scale,'LineWidth',2);
datetick('x','yy','keeplimits');
grid on; axis tight;
title('CIP','interpreter','latex','Fontsize',15);
subplot(3,3,8);plot(dates(datesperiod),(ois(datesperiod))*abs_scale,'LineWidth',2);
datetick('x','yy','keeplimits');
grid on; axis tight;
title('OIS','interpreter','latex','Fontsize',15);
subplot(3,3,9);plot(dates(datesperiod),(Chi_D_US(datesperiod))*abs_scale,'LineWidth',2);
datetick('x','yy','keeplimits');
grid on; axis tight;
title('$\chi^-_{us}-\chi^+_{us}$','interpreter','latex','Fontsize',15);
h=suptitle('Data');
set(h,'interpreter','latex','Fontsize',20);
if printit==1
    orient landscape;
    set(gcf, 'PaperUnits', 'normalized');
    set(gcf, 'PaperPosition', [0 0 1 1]);  
    print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
    print(gcf,'-dpdf',[foldername 'raw_data']);
end

%% Figure w-6 Panels
figure
subplot(3,2,1);
plot(dates(datesperiod),exp(-inv_e(datesperiod)),'LineWidth',2);
datetick('x','yy','keeplimits');
grid on; axis tight;
title('$e $ (EUR/USD)','interpreter','latex','Fontsize',15);
subplot(3,2,2);
plot(dates(datesperiod),exp(mu_us(datesperiod)),'LineWidth',2); hold on;
plot(dates(datesperiod),exp(mu_eu(datesperiod)),'LineWidth',2); hold off;
h=legend('c=US','c=EU','color','none','box','off','location','northwest');
set(h,'interpreter','latex');
datetick('x','yy','keeplimits');
grid on; axis tight;
title('$\mu^{c}$','interpreter','latex','Fontsize',15);
subplot(3,2,3);
plot(dates(datesperiod),exp(im_eu(datesperiod)),'LineWidth',2); hold on;
plot(dates(datesperiod),exp(im_us(datesperiod)),'LineWidth',2);
datetick('x','yy','keeplimits'); 
grid on; axis tight;
title('$i_m^{c}$','interpreter','latex','Fontsize',15);
h=legend('c=US','c=EU','color','none','box','off','location','northeast');
set(h,'interpreter','latex');
subplot(3,2,4);
plot(dates(datesperiod),exp(M_us(datesperiod)),'LineWidth',2); hold on;
plot(dates(datesperiod),exp(M_eu(datesperiod)),'LineWidth',2); hold off;
h=legend('c=US','c=EU','color','none','box','off','location','northwest');
set(h,'interpreter','latex');
datetick('x','yy','keeplimits');
grid on; axis tight;
title('$M^{c}$','interpreter','latex','Fontsize',15);
subplot(3,2,5);
plot(dates(datesperiod),Rb_Rm(datesperiod)*abs_scale,'LineWidth',2); hold on;
plot(dates(datesperiod),cip(datesperiod)*abs_scale,'LineWidth',2); 
plot(dates(datesperiod),Rb_Rm_eu(datesperiod)*abs_scale,'LineWidth',2); hold off;
h=legend('$\mathcal{BP}$','$\mathcal{CIP}$','$\mathcal{BP_{eu}}$','color','none','box','off','location','northwest');
set(h,'interpreter','latex');
datetick('x','yy','keeplimits');
grid on; axis tight;
%title('$R_b^*-R_m^*$','interpreter','latex','Fontsize',15);
%title('CIP','interpreter','latex','Fontsize',15);
title('$\mathcal{BP},\mathcal{CP}$','interpreter','latex','Fontsize',15);
subplot(3,2,6);
plot(dates(datesperiod),ois(datesperiod)*abs_scale,'LineWidth',2); hold on;
plot(dates(datesperiod),Chi_D_US(datesperiod)*abs_scale,'LineWidth',2);
datetick('x','yy','keeplimits');
grid on; axis tight;
title('Other Spreads','interpreter','latex','Fontsize',15);
h=legend('$\mathcal{OIS}$','$R^f Spread $','color','none','box','off','location','northwest');
set(h,'interpreter','latex');
h=suptitle('Data');
set(h,'interpreter','latex','Fontsize',20);
if printit==1
    % orient landscape;
    set(gcf, 'PaperUnits', 'normalized');
    set(gcf, 'PaperPosition', [0 0 1 1]);  
  %  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
    print(gcf,'-dpdf',[foldername 'raw_data']);
end