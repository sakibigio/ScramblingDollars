%% Liquidity Exchange Rate Model: Get Data and Predetermined Estimations
% This file reads data from raw files, calculates steady-state values of
% observables, and estimate the AR sequences of interest rates and money
% bases outside of Dynare
clear; 
close all;

% Notes: July '21
% issues with persistence...we need to externally calibrate persistence.

%% List of model variables
% Log of Liquidity ratio: log_liqratio mu_us mu_eu
% CIP, LIBOR and OIS data: ois cip libor gap
% Policy Rates

% Target Years
year_target=(25:72);
exo_persistence_M=1; % Exogenous persistence of M
exo_persistence=0  ; % Exogenous persistence of i

% Inflation and Monetary Base
plotit=0;
dates=datenum(2001,1:234,1);
dates=datenum(2001,1:234,1);
formataxis   = @(x) set(x, 'Fontname', 'Times', 'FontWeight', 'normal');
label_x      = @(x) xlabel(x,'Fontname', 'Times', 'FontWeight', 'normal','Interpreter','latex');
label_y      = @(x) ylabel(x,'Fontname', 'Times', 'FontWeight', 'normal','interpreter','latex');
    
% Adjustments
freq = 12;
liq_ratio_scale=0.18     ; % Scale used to suggest ratio
liq_ratio_eu_scale=0.40  ; % Scale used to match a bank ratio

% monthly decimal return to annual BPS
abps_factor=freq*1e4;

% Detrend Money Supply
detrendM=1;
demean_range=year_target;

%% -> we add a scaling factor to the EBP, else we cannot get Rb-Rm*>Rm-Rm*...
Rb_Rm_scale=0/freq/1e4; % Shift in the scale of Rb-Rm -> like risk premium
% Rb_Rm_scale=200/freq/1e4; % Shift in the scale of Rb-Rm -> like risk premium
iota_ss = 0.1;            % Discount window-IOR spread (decimal, annual, nominal)


%% US liquidity ratio
log_liqratio_t = log(liq_ratio_scale*xlsread('LFX_datainputs.xlsx','DataCounterpart','C27:C260'));

% % Option: we can detrend liquidity ratio
% ln_mu_us_obs = diff(log_liqratio_t);
% ln_mu_us_obs = ln_mu_us_obs-mean(ln_mu_us_obs);
display(['Average liquidity ratio ' num2str(mean(exp(log_liqratio_t)))]);

% Currently we do not detrend it
%mu_us = log_liqratio_t(1)+cumsum([0;ln_mu_us_obs]);
mu_us = log_liqratio_t;

%% Euro Liquidity Ratio
temp  = xlsread('LFX_datainputs.xlsx','DataCounterpart','D27:D260');
mu_eu = log(liq_ratio_eu_scale*temp);

%% Spreads data: OIS, Libor Gap, CIP
cip      = xlsread('LFX_datainputs.xlsx','DataCounterpart','E27:E260')/freq/1e4;
ois      = xlsread('LFX_datainputs.xlsx','DataCounterpart','F27:F260')/freq/1e4;
ebp      = xlsread('LFX_datainputs.xlsx','DataCounterpart','G27:G260')/freq/1e4;
ebp_eu   = xlsread('LFX_datainputs.xlsx','DataCounterpart','CG27:CG260')/freq/1e4;
Chi_D_US = xlsread('LFX_datainputs.xlsx','DataCounterpart','H27:H260')/freq/1e4;


%% Policy interest rates
% Policy Rates
policyrate = xlsread('LFX_datainputs.xlsx','DataCounterpart','S27:AB260');
policyrate(isnan(policyrate)) = 0; % Check issues here
Liborrate = xlsread('LFX_datainputs.xlsx','DataCounterpart','AI27:AQ260');
Liborrate(isnan(Liborrate)) = 0; % Check issues here
TED_all = xlsread('LFX_datainputs.xlsx','DataCounterpart','CQ27:CZ260');
TED_s_us_t = TED_all(:,1)/freq/100;
i_us_t = policyrate(:,1)/freq/100;
BP_s_all = xlsread('LFX_datainputs.xlsx','DataCounterpart','CG27:CP260');
BP_s_us_t = BP_s_all(:,1)/freq/100;
CIP_s_all = xlsread('LFX_datainputs.xlsx','DataCounterpart','DA27:DJ260');
CIP_s_us_t = CIP_s_all(:,1)/freq/100;

%% Discount Window Counterpart
DW_t = xlsread('LFX_datainputs.xlsx','DataCounterpart','BV27:BV260'); % Only Primary Credit

%% Inflation and Money Base
% US money base M_us
%res_us_t = xlsread('LFX_rawdata.xls','EuroData','T28:T261');
%sec_us_t = xlsread('LFX_rawdata.xls','EuroData','U28:U261');
M_us_t = xlsread('LFX_datainputs.xlsx','DataCounterpart','AC27:AC260');
liq_t = (M_us_t);
log_liq_t = log(liq_t);
%log_liq_t_detrend = log_liq_t-(1:length(log_liq_t))'*log(mean(1+pi_us_t));
d_log_liq_t = diff(log_liq_t);

% Detrend Money Options
if detrendM==1
    gM_us       = mean(d_log_liq_t(demean_range));
    d_log_liq_t = d_log_liq_t-gM_us;
    M_us = log_liq_t(1)+cumsum([0; d_log_liq_t]);
    M_us_ss = exp(mean(M_us));
    X = [ones(length(M_us)-1,1) M_us(1:end-1)];
    [B,~,R,~,~] = regress(M_us(2:end),X);
    rho_M_us = B(2);
    sigma_M_us = std(R);
else
    d_log_liq_t = d_log_liq_t;
    M_us = log_liq_t(1)+cumsum([0; d_log_liq_t]);
    M_us_ss = exp(mean(M_us));
    X = [ones(length(M_us)-1,1) M_us(1:end-1)];
    [B,~,R,~,~] = regress(M_us(2:end),X);
    rho_M_us = B(2);
    sigma_M_us = std(R);
end
    
% Definition of M
M_us = log_liq_t(1)+cumsum([0; d_log_liq_t]);
M_us_ss = exp(mean(M_us));
X = [ones(length(M_us)-1,1) M_us(1:end-1)];
[B,~,R,~,~] = regress(M_us(2:end),X);
rho_M_us = B(2);
sigma_M_us = std(R);

% Euro Money
M_eu_t         = xlsread('LFX_datainputs.xlsx','DataCounterpart','AD28:AD261');
liq_eu_t       = M_eu_t; 
log_liq_eu_t   = log(M_eu_t);
d_log_liq_eu_t = diff(log_liq_eu_t);
if detrendM==1   
    gM_eu       = mean(d_log_liq_eu_t(demean_range));
    d_log_liq_t = d_log_liq_t-gM_eu;
    d_log_liq_eu_t = d_log_liq_eu_t-mean(d_log_liq_eu_t);
    M_eu = log_liq_eu_t(1)+cumsum([0;d_log_liq_eu_t]);
    M_eu_ss = exp(mean(M_eu));
    X = [ones(length(M_eu)-1,1) M_eu(1:end-1)];
    [B,~,R,~,~] = regress(M_eu(2:end),X);
    rho_M_eu = B(2);
    sigma_M_eu = std(R);
else
    d_log_liq_eu_t = d_log_liq_eu_t;
    M_eu = log_liq_eu_t(1)+cumsum([0;d_log_liq_eu_t]);
    M_eu_ss = exp(mean(M_eu));
    X = [ones(length(M_eu)-1,1) M_eu(1:end-1)];
    [B,~,R,~,~] = regress(M_eu(2:end),X);
    rho_M_eu = B(2);
    sigma_M_eu = std(R);    
end    

if plotit==1
    figure;
    plot(dates,M_us_t); hold on; plot(dates,M_eu_t); hold off; pause(1)
    legend('US','EU');
    datetick('x','yy');
    grid on; axis tight;
    label_x('Time (Year)');
    title('Money Supply');
end
% persistence of Euro Money Supply
if exo_persistence_M==1
    rho_M_us=0.995;
    rho_M_eu=0.995;
end


%% Other Currencies
% US inflation rate pi_us and policy interest rate im_us
pi_us_t = xlsread('LFX_datainputs.xlsx','EuroData','AC28:AC261')/freq/100;
libor_us_t = xlsread('LFX_datainputs.xlsx','EuroData','K28:K261')/freq/100;

% Other countries' policy rate im_co, inflation rate pi_co and exchange
% rate against USD
curlist = {'au','ca','eu','jp','nz','no','sw','ch','uk','us'};
datalist = {'Australia','Canada','Euro','Japan','NewZealand','Norway','Sweden','SwissFranc','UK'};
curindexlist = [3,2,10,9,4,6,8,7,5];

% Read data of policy rates, inflation and exchange rates
for j=1:length(datalist)
%     eval(['i_' curlist{j} '_t = xlsread(''LFX_rawdata.xls'',''' datalist{j} 'Data'',''L28:L261'')/freq/100;']);
   % eval(['libor_' curlist{j} '_t = xlsread(''LFX_datainputs.xlsx'',''' datalist{j} 'Data'',''L28:L261'')/freq/100;']);
    eval(['i_' curlist{j} '_t     = policyrate(:,' num2str(curindexlist(j)) ')/freq/100;']);
    eval(['libor_' curlist{j} '_t     = policyrate(:,' num2str(curindexlist(j)) ')/freq/100;']);
    eval(['TED_s_' curlist{j} '_t     = TED_all(:,' num2str(curindexlist(j)) ')/freq/100;']);
    eval(['pi_' curlist{j} '_t    = xlsread(''LFX_datainputs.xlsx'',''' datalist{j} 'Data'',''AD28:AD261'')/freq/100;']);
    eval(['ln_' curlist{j} '_us_t = xlsread(''LFX_datainputs.xlsx'',''' datalist{j} 'Data'',''C28:C261'');']);
    eval(['BP_s_' curlist{j} '_t     = BP_s_all(:,' num2str(curindexlist(j)) ')/freq/100;']);
    eval(['CIP_s_' curlist{j} '_t     = CIP_s_all(:,' num2str(curindexlist(j)) ')/freq/100;']);
end

% Calculate steady state policy rates and AR(1) parameters
savelist = '';
for j=1:length(curlist)
    eval(['i_t_temp = i_' curlist{j} '_t;']);
    eval(['libor_t_temp = libor_' curlist{j} '_t;']);
    eval(['imss_' curlist{j} '=mean(1+i_t_temp)/mean(1+pi_' curlist{j} '_t);']);
    eval(['iwss_' curlist{j} '=mean(1+i_t_temp+iota_ss/freq)/mean(1+pi_' curlist{j} '_t);']);
    eval(['RLiborss_' curlist{j} '=mean(1+libor_t_temp)/mean(1+pi_' curlist{j} '_t);']);
    eval(['i_t_temp = (1+i_t_temp)/mean(1+pi_' curlist{j} '_t)-1;']);
    eval(['libor_t_temp = (1+libor_t_temp)/mean(1+pi_' curlist{j} '_t)-1;']);
    X=[ones(length(i_t_temp)-1,1) i_t_temp(1:end-1)];
    [B,BINT,R,RINT,STATS] = regress(i_t_temp(2:end),X);
    eval(['i_' curlist{j} '_t = i_t_temp+1;']);
    eval(['libor_' curlist{j} '_t = libor_t_temp+1;']);
    %% -> CHECKOUT PERSISTENCE
    eval(['rho_im_' curlist{j} '=B(2);']);
    if exo_persistence==1
        eval(['rho_im_' curlist{j} '=0.99;']);
    end
    eval(['sigma_im_' curlist{j} '=std(R);']);
    eval(['piss_' curlist{j} '=mean(pi_' curlist{j} '_t)+1;']);
    savelist = [savelist 'imss_' curlist{j} ' iwss_' curlist{j} ' rho_im_' curlist{j} ' sigma_im_' curlist{j} ' RLiborss_' curlist{j} ' '];
end

% Subtract 100bps from AUS policy rate and 200bps from NZL policy rates
imss_au = mean(i_au_t)-100/12/1e4;
imss_nz = mean(i_nz_t)-200/12/1e4;
imss_no = mean(i_no_t)-50/12/1e4;
iwss_au = iwss_au-100/12/1e4;
iwss_nz = iwss_nz-200/12/1e4;
iwss_no = iwss_no-50/12/1e4;
RLiborss_au = mean(libor_au_t)-100/12/1e4;
RLiborss_nz = mean(libor_nz_t)-200/12/1e4;
RLiborss_no = mean(libor_no_t)-50/12/1e4;

% only affects rates

% Save the estimated parameters of money base and policy rates for Dynare
% -> saves US and Euro Money Supply

%% -> CHECKOUT PERSISTENCE OF RATES!
% eval(['rho_im_' curlist{j} '=0.991;']);
% rho_im_us=rho_im_jp-0.007;
% rho_eu_us=rho_im_eu+0.002;
% rho_im_us=0.999;
% rho_im_eu=0.999;
eval(['save dynare_calibration_param.mat M_us_ss rho_M_us sigma_M_us M_eu_ss rho_M_eu sigma_M_eu ' savelist ';']);

%% Saving FX Data
% Calculate and save exchange rates (steady states and time-series paths)
savelist = '';
for j=1:length(curlist)-1
%    eval(['ln_' curlist{j} '_us_t = ln_' curlist{j} '_us_t-(1:length(ln_' curlist{j} '_us_t))''*log(mean(1+pi_' curlist{j} '_t)/mean(1+pi_us_t));']);
    eval(['ln_' curlist{j} '_us_ss = mean(ln_' curlist{j} '_us_t);']);
    eval(['inv_e_' curlist{j} ' = -ln_' curlist{j} '_us_t;']);
    savelist = [savelist 'ln_' curlist{j} '_us_t ln_' curlist{j} '_us_ss '];
end
inv_e = inv_e_eu;
eval(['save exchange_rate_data.mat ' savelist ';']);

% Save the time series of policy rates
%savelist = '';
for j=1:length(curlist)
    eval(['im_' curlist{j} ' = log(i_' curlist{j} '_t);']);
    eval(['RLibor_' curlist{j} ' = log(libor_' curlist{j} '_t);']);
    eval(['im_' curlist{j} '_obs = diff(log(i_' curlist{j} '_t));']);
    savelist = [savelist 'im_' curlist{j} '_obs im_' curlist{j} ' RLibor_' curlist{j} ' TED_s_' curlist{j} '_t ' ' BP_s_' curlist{j} '_t ' ...
        ' CIP_s_' curlist{j} '_t '];
end

%% -> Subtract 100bps from AUS policy rate and 200bps from NZL policy rates
im_au = im_au-100/12/1e4;
im_nz = im_nz-200/12/1e4;
RLibor_au = RLibor_au-100/12/1e4;
RLibor_nz = RLibor_nz-200/12/1e4;

%% Excess Bond Premium
Rb_Rmtemp = ebp;
Rb_Rmtemp_eu = ebp_eu;
Rb_Rm     = Rb_Rmtemp-mean(Rb_Rmtemp)*(Rb_Rm_scale>0)+Rb_Rm_scale         ; % There's an adjustment here...
Rb_Rm_eu  = Rb_Rmtemp_eu-mean(Rb_Rmtemp_eu)*(Rb_Rm_scale>0)+Rb_Rm_scale; % There's an adjustment here...
Rb_us     = log(Rb_Rm+exp(im_us));

% if plotit==1
%     figure;
%     plot(cip); hold on; plot(ois); plot(ebp); plot(Chi_D_US); hold off;
%     legend('Bond CIP','OIS CIP','Libor','Chi Diff.');
% end

%% Pre Great Recession Data
Ted_us_yt=mean(TED_s_us_t(year_target))*1e4*12;
Ted_eu_yt=mean(TED_s_eu_t(year_target))*1e4*12;
ois_yt=mean(ois(year_target))*1e4*12;
cip_yt=mean(cip(year_target))*1e4*12;
uip_yt=mean(i_eu_t(year_target)-i_us_t(year_target))*1e4*12;
mu_us_yt=mean(exp(mu_us(year_target)));
mu_eu_yt=mean(exp(mu_eu(year_target)));
Rb_Rm_yt=mean(Rb_Rm(year_target)*freq*1e4);
inv_e_yt=mean(exp(inv_e(year_target)));
save LFX_targets Ted_us_yt Ted_eu_yt mu_us_yt mu_eu_yt Rb_Rm_yt inv_e_yt ois_yt cip_yt uip_yt;

%% Conversions to desired format
bps_scale=12e4;
ois_ts   = ois*bps_scale;
cip_ts   = cip*bps_scale;
idiff_ts = (i_eu_t-i_us_t)*bps_scale;
Rb_Rm_ts = Rb_Rm*bps_scale;
mu_us_ts = exp(mu_us);
mu_eu_ts = exp(mu_eu);
e_euus_ts= exp(inv_e);
% Variable
var_list={'ois','cip','idiff','Rb_Rm','mu_us','mu_eu','e_euus'};
mom_list=[];
for j=1:numel(var_list)
     eval(['aux=' var_list{j} '_ts;']);
     rho= autocorr(aux);
     m  = mean(aux(year_target)); %subsample
     s  = std(aux); 
     eval([var_list{j} '_m.mean=' num2str(m) ';']);
     eval([var_list{j} '_m.std=' num2str(s)  ';']);
     eval([var_list{j} '_m.rho=' num2str(rho(2)) ';']);
     mom_list = [mom_list ' ' var_list{j} '_m'];
end
eval(['save LFX_datamoments.mat' mom_list]);
clear aux;
% moment_list={'mean','std','autocorr'}
% for 
% for j=1:length(curlist)
% 
% end

%% Some Plots
% Plot Deviations
if plotit==1
    aux=im_eu-im_us;
    figure;
    plot(dates,cip*abps_factor); hold on; 
    plot(dates,ois*abps_factor); 
    plot(dates,ebp*abps_factor); 
    plot(dates,aux*abps_factor); 
    legend('Bond CIP','OIS CIP','EBP','interest differential');
 %   legend('Bond CIP','OIS CIP','EBP','US Differential');
%     plot(dates,Chi_D_US*abps_factor); hold off;
%    legend('Bond CIP','OIS CIP','EBP','Chi Diff.');
    datetick('x','yy');
    grid on; axis tight;
    label_x('Time (Year)');
end

if plotit==1
    figure;
    plot(dates,Rb_Rm); hold on; 
    plot(dates,Rb_Rm_eu); 
    legend('Bond Premium','Euro Bond Premium');
 %   legend('Bond CIP','OIS CIP','EBP','US Differential');
%     plot(dates,Chi_D_US*abps_factor); hold off;
%    legend('Bond CIP','OIS CIP','EBP','Chi Diff.');
    datetick('x','yy');
    grid on; axis tight;
    label_x('Time (Year)');
end

if plotit==1
    figure
    plot(dates,exp(mu_us)); hold on; plot(dates,exp(mu_eu)); hold off;  grid on;
    legend('US','EU');
    datetick('x','yy');
    grid on; axis tight;
    label_x('Time (Year)'); pause(1);
end

if plotit==1
    figure;
    plot(dates,exp(M_us)/M_us_ss); hold on; plot(dates,exp(M_eu)/M_eu_ss); hold off; 
    legend('US','EU');
    datetick('x','yy');
    grid on; axis tight;
    label_x('Time (Year)');
    title('Detrended Money Supply'); pause(1)

    figure;
    plot(dates,exp(M_us)/M_us_ss); hold on; plot(dates,exp(mu_us)/mu_us_yt); hold off; 
    legend('M US','Liquidity Ratio');
    datetick('x','yy');
    grid on; axis tight;
    label_x('Time (Year)');
    title('Detrended Money Supply vs. Liquidity Ratio'); pause(1);
    
    figure;
    plot(dates,exp(M_eu)/M_eu_ss); hold on; plot(dates,exp(mu_eu)/mean(exp(mu_eu))); hold off; 
    legend('M EU','Liquidity Ratio');
    datetick('x','yy');
    grid on; axis tight;
    label_x('Time (Year)');
    title('Detrended Money Supply vs. Liquidity Ratio'); pause(1);

    figure
    plot(dates,(exp(M_us)/M_us_ss)./(exp(M_eu)/M_eu_ss)); 
    datetick('x','yy');
    grid on; axis tight;
    label_x('Time (Year)');
    title('Ratio of Detrended Money Supply'); pause(1);
end

%% Save for Estimation
Ted_us=TED_s_us_t;
Ted_eu=TED_s_eu_t;

%% Save all data for dynare estimation
eval(['save LFX_data3.mat mu_eu mu_us inv_e Ted_us Ted_eu inv_e_jp inv_e_ch M_us Rb_Rm Rb_Rm_eu Rb_us M_eu ois cip Chi_D_US pi_us_t pi_eu_t DW_t ' savelist ';']);
