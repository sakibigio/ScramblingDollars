% compute_moments.m
% Extracted from main_LFX.m (lines 2491–2695) on 2026-02-28
% Called by: main_LFX.m
%
% Contents:
%   - Theoretical moments (means, std devs by regime)
%   - Regime-switching deviations
%   - Simulation-based moments (means, stds, autocorrelations)
%   - LaTeX table output: Mod_SwitchDev_Moments.tex, 
%     Mod_Moments.tex, Mod_CIP_Moments.tex

%% Moments - Theoretical
bp_vec=(Rb_us_vec.^(freq)-Rm_us_vec.^(freq))*rate_scale;
cip_vec=(Rm_eu_vec.^(freq)-Rm_us_vec.^(freq))*rate_scale;
pi_us_ret_vec=(pi_us_vec.^(freq)-1)*rate_scale;

% CIP check
figure('Name','CIP+Ret comp','NumberTitle','off');
plot(sigma_us_vec(index1),0*Rm_eu_vec(index1), 'LineStyle', '--','Color','k','LineWidth',1); hold on;
splot1(sigma_us_vec(index1),cip_vec(index1)); hold on;
splot2(sigma_us_vec(index1),pi_us_ret_vec(index1));
splot3(sigma_us_vec(index2),cip_vec(index2)); 
splot4(sigma_us_vec(index2),pi_us_ret_vec(index2));
hold off;
%splot3(sigma_us_vec(index2),diff);
xlim([sigma_us_vec(1) sigma_us_vec(end)]);

% Construct the regime specific vectors
e_euus_vec_r1=e_euus_vec(index1);
bp_vec_r1=bp_vec(index1);
cip_vec_r1=cip_vec(index1);
pi_us_ret_vec_r1=pi_us_ret_vec(index1);
e_euus_vec_r2=e_euus_vec(index2);
bp_vec_r2=bp_vec(index2);
cip_vec_r2=cip_vec(index2);
pi_us_ret_vec_r2=pi_us_ret_vec(index2);

% means:
E_e=e_euus_vec*invp;
E_bp=bp_vec*invp;
E_cip=cip_vec*invp;
E_pi_ret_us=pi_us_ret_vec*invp;
E_e_r1=e_euus_vec_r1*invp1;
E_bp_r1=bp_vec_r1*invp1;
E_cip_r1=cip_vec_r1*invp1;
E_pi_ret_us_r1=pi_us_ret_vec_r1*invp1;
E_e_r2=e_euus_vec_r2*invp2;
E_bp_r2=bp_vec_r2*invp2;
E_cip_r2=cip_vec_r2*invp2;
E_pi_ret_us_r2=pi_us_ret_vec_r2*invp2;

% standard deviations:
std_e=sqrt((e_euus_vec-E_e).^2*invp);
std_bp=sqrt((bp_vec-E_bp).^2*invp);
std_cip=sqrt((cip_vec-E_cip).^2*invp);
std_pi_ret_us=sqrt((pi_us_ret_vec-E_pi_ret_us).^2*invp);
std_e_r1=sqrt((e_euus_vec_r1-E_e_r1).^2*invp1);
std_bp_r1=sqrt((bp_vec_r1-E_bp_r1).^2*invp1);
std_cip_r1=sqrt((cip_vec_r1-E_cip_r1).^2*invp1);
std_pi_ret_us_r1=sqrt((pi_us_ret_vec_r1-E_pi_ret_us_r1).^2*invp1);
std_e_r2=sqrt((e_euus_vec_r2-E_e_r2).^2*invp2);
std_bp_r2=sqrt((bp_vec_r2-E_bp_r2).^2*invp2);
std_cip_r2=sqrt((cip_vec_r2-E_cip_r2).^2*invp2);
std_pi_ret_us_r2=sqrt((pi_us_ret_vec_r2-E_pi_ret_us_r2).^2*invp2);

% Moments Conditional on Jump
Dev_Mat=(((e_euus_vec').^(-1)*e_euus_vec)-1)*rate_scale;
Jump_mat=(Zprob_sigma_us.*Dev_Mat);
Jump_test=Zprob_sigma_us;
Dev_r1_to_r2=invp1'*(Jump_mat(1:N_sigma_us/2,N_sigma_us/2+1:end))/P(1,2)*ones(N_sigma_us/2,1);
%Dev_test=invp1'*(Jump_test(1:N_sigma_us/2,N_sigma_us/2+1:end))/P(1,2)*ones(N_sigma_us/2,1);
Dev_r2_to_r1=invp2'*(Jump_mat(N_sigma_us/2+1:end,1:N_sigma_us/2))/P(2,1)*ones(N_sigma_us/2,1);

filename = fullfile(foldername, 'Mod_SwitchDev_Moments.tex');
fid = fopen(filename, 'wt');
% devswith_moments = {...
    %'FX', 1, rho_e_sim, std_e/E_e, (E_e_r2./E_e_r1-1)*10000, rho_e_sim_r2, rho_e_sim_r1, (std_e_r2./std_e_r1); 
    %             '$\Delta$ FX', E_pi_ret_us, rho_pi_us_sim, std_pi_ret_us, E_pi_ret_us_r2-E_pi_ret_us_r1, rho_e_sim_r2, rho_e_sim_r1, (std_pi_ret_us_r2./std_pi_ret_us_r1);
    %             'BP', E_bp, rho_bp_sim, std_bp, E_bp_r2-E_bp_r1, rho_bp_sim_r2, rho_bp_sim_r1, (std_bp_r2./std_bp_r1);                 
               % , };

% fprintf(fid, '\\bottomrule\n');
%for ii = 1:size(model_moments, 1)
    fprintf(fid, '%.1f & %.1f \\\\', Dev_r1_to_r2,Dev_r2_to_r1);
%end
fclose(fid);

%% Moments - Simulation Based
% means:
E_e_sim=mean(e_euus_t);
E_bp_sim=mean(bp_us_t);
E_cip_sim=mean(uip_Rm_t);

% standard deviations:
std_e_sim=std(e_euus_t);
std_bp_sim=std(bp_us_t);
std_cip_sim=std(uip_Rm_t);

% autocorrelation:
aux=autocorr(e_euus_t);
rho_e_sim=aux(2);
aux=autocorr(bp_us_t);
rho_bp_sim=aux(2);
aux=autocorr(uip_Rm_t);
rho_cip_sim=aux(2);
aux=autocorr(pi_us_t);
rho_pi_us_sim=aux(2);

% Autocorrelation by Regime
index_r1=chain<=N_sigma_us/2;
index_r2=chain>N_sigma_us/2;

% e, bp, cip
cip_t=uip_Rm_t;
aux=autocorr(e_euus_t(index_r1));
rho_e_sim_r1=aux(2);
aux=autocorr(e_euus_t(index_r2));
rho_e_sim_r2=aux(2);
aux=autocorr(bp_us_t(index_r1));
rho_bp_sim_r1=aux(2);
aux=autocorr(bp_us_t(index_r2));
rho_bp_sim_r2=aux(2);
aux=autocorr(cip_t(index_r1));
rho_cip_sim_r1=aux(2);
aux=autocorr(cip_t(index_r2));
rho_cip_sim_r2=aux(2);
aux=autocorr(pi_us_t(index_r1));
rho_pi_us_sim_r1=aux(2);
aux=autocorr(pi_us_t(index_r2));
rho_pi_us_sim_r2=aux(2);



%% Print Model Moments to Table
% Define folder and file name
% filename = fullfile(foldername, 'Mod_Moments.tex');
% fid = fopen(filename, 'wt');

% Write the LaTeX table header
% fprintf(fid, '\\begin{table}[h!]\n');
% fprintf(fid, '\\centering\n');
% fprintf(fid, '\\begin{threeparttable}\n');
% fprintf(fid, '\\begin{tabular}{lcccccc}\n');
% fprintf(fid, '\\toprule\n');
% fprintf(fid, '\\textbf{Variable} & \\multicolumn{6}{c}{\\textbf{Moment}} \\\\\n');
% fprintf(fid, '\\cmidrule(lr){2-7}\n');
% fprintf(fid, ' & \\multicolumn{3}{c}{Unconditional} & \\multicolumn{3}{c}{Regime 1 / Regime 2} \\\\\n');
% fprintf(fid, '\\cmidrule(lr){2-4} \\cmidrule(lr){5-7}\n');
% fprintf(fid, ' & Mean & Autocorrelation & Std & Mean & Autocorrelation & Std \\\\\n');
% fprintf(fid, '\\midrule\n');

% % Section for Data Moments
% fprintf(fid, '\\rowcolor{gray!20}\n');
% fprintf(fid, '\\multicolumn{7}{c}{\\textbf{Data Moments}} \\\\\n');

% % Replace the following with your actual data
% data_moments = {'FX', 0.123, 0.456, 0.789, 0.123, 0.456, 0.789; 
%                 'BP', 0.223, 0.556, 0.889, 0.223, 0.556, 0.889;
%                 'CIP', 0.323, 0.656, 0.989, 0.323, 0.656, 0.989};

% for ii = 1:size(data_moments, 1)
%     fprintf(fid, '%s & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f \\\\\n', data_moments{ii, :});
% end

% Section for Model Moments
% fprintf(fid, '\\midrule\n');
% fprintf(fid, '\\rowcolor{gray!20}\n');
% fprintf(fid, '\\multicolumn{7}{c}{\\textbf{Model Moments}} \\\\\n');

% Replace the following with your actual model moments
filename = fullfile(foldername, 'Mod_Moments.tex');
fid = fopen(filename, 'wt');
model_moments = {'FX', 1, rho_e_sim, std_e/E_e, (E_e_r2./E_e_r1-1)*10000, rho_e_sim_r2, rho_e_sim_r1, (std_e_r2./std_e_r1); 
                 '$\Delta$ FX', E_pi_ret_us, rho_pi_us_sim, std_pi_ret_us, E_pi_ret_us_r2-E_pi_ret_us_r1, rho_e_sim_r2, rho_e_sim_r1, (std_pi_ret_us_r2./std_pi_ret_us_r1);
                 'BP', E_bp, rho_bp_sim, std_bp, E_bp_r2-E_bp_r1, rho_bp_sim_r2, rho_bp_sim_r1, (std_bp_r2./std_bp_r1);                 
                 'CIP', E_cip, rho_cip_sim, std_cip, E_cip_r2-E_cip_r1, rho_cip_sim_r2, rho_cip_sim_r1, (std_cip_r2./std_cip_r1)};

% fprintf(fid, '\\bottomrule\n');
for ii = 1:size(model_moments, 1)
    fprintf(fid, '%s & %.1f & %.2f & %.1f & %.1f &  \\{ %.2f, %.2f \\} & %.1f \\\\ \n', model_moments{ii, :});
end

% Close the table
% fprintf(fid, '\\bottomrule \n');
fclose(fid);

filename = fullfile(foldername, 'Mod_CIP_Moments.tex');
fid = fopen(filename, 'wt');
model_moments = {...
    %'FX', 1, rho_e_sim, std_e/E_e, (E_e_r2./E_e_r1-1)*10000, rho_e_sim_r2, rho_e_sim_r1, (std_e_r2./std_e_r1); 
    %             '$\Delta$ FX', E_pi_ret_us, rho_pi_us_sim, std_pi_ret_us, E_pi_ret_us_r2-E_pi_ret_us_r1, rho_e_sim_r2, rho_e_sim_r1, (std_pi_ret_us_r2./std_pi_ret_us_r1);
    %             'BP', E_bp, rho_bp_sim, std_bp, E_bp_r2-E_bp_r1, rho_bp_sim_r2, rho_bp_sim_r1, (std_bp_r2./std_bp_r1);                 
                 'CIP (model)', E_cip, rho_cip_sim, std_cip, E_cip_r2-E_cip_r1, rho_cip_sim_r2, rho_cip_sim_r1, (std_cip_r2./std_cip_r1)};

% fprintf(fid, '\\bottomrule\n');
for ii = 1:size(model_moments, 1)
    fprintf(fid, '%s & %.1f & %.2f & %.1f & %.1f &  \\{ %.2f, %.2f \\} & %.1f \n', model_moments{ii, :});
end

% Close the table
% fprintf(fid, '\\bottomrule \n');
fclose(fid);
% fprintf(fid, '\\end{tabular}\n');

% % Add the units and notes
% fprintf(fid, '\\begin{tablenotes}\n');
% fprintf(fid, '\\small\n');
% fprintf(fid, '\\item \\textbf{Units:} \\newline Mean: Percentage (\\%%), Autocorrelation: Unitless, Std: Percentage (\\%%)\n');
% fprintf(fid, '\\end{tablenotes}\n');
% fprintf(fid, '\\end{threeparttable}\n');
% fprintf(fid, '\\caption{Data vs. Model Moments.}\n');
% fprintf(fid, '\\label{tab:example}\n');
% fprintf(fid, '\\end{table}\n');

% Close the file

