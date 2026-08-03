%% plot_expected_inflation.m
% Figure for the G07 verification objects computed by run_expected_inflation.m.
% Loads data/expected_inflation_cd.mat (grids over the sigma_us Markov chain).
%
% Panel (a): expected inflation E[pi'|s], annualized %, by currency and regime,
%            with the ergodic distribution of sigma_us shaded underneath.
% Panel (b): Jensen wedge  E[Pi'] - 1/E[1/Pi']  in monthly bps, state by state
%            (the convention gap flagged in audit item G07).
%
% Invoke:  matlab -batch "addpath('plotting'); plot_expected_inflation"

S = load(fullfile('data','expected_inflation_cd.mat'));
N_s  = numel(S.Epi_us_grid);
i_r1 = 1:N_s/2;  i_r2 = N_s/2+1:N_s;
sig  = S.sigma_us_vec(i_r1);            % same sigma nodes in both regimes
ann  = @(g) (g.^S.freq - 1)*100;
bpsm = @(g) g*S.rate_scale;

cUS = [0.85 0.2 0.2]; cEU = [0.1 0.35 0.7];

fig = figure('Position',[100 100 1150 430],'Color','w');

% ---------- (a) expected inflation ----------
subplot(1,2,1);
yyaxis right;                            % ergodic density in the background
erg1 = S.erg(i_r1); erg2 = S.erg(i_r2);
hE1 = area(sig, erg1,'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.25,'EdgeColor','none'); hold on;
hE2 = area(sig, erg2,'FaceColor',[0.9 0.6 0.2],'FaceAlpha',0.25,'EdgeColor','none');
ylabel('Ergodic probability','Interpreter','latex');
set(gca,'YColor',[0.4 0.4 0.4]);
yyaxis left;
hU1 = plot(sig, ann(S.Epi_us_grid(i_r1)),'-','Color',cUS,'LineWidth',2.2); hold on;
hU2 = plot(sig, ann(S.Epi_us_grid(i_r2)),'--','Color',cUS,'LineWidth',2.2);
hV1 = plot(sig, ann(S.Epi_eu_grid(i_r1)),'-','Color',cEU,'LineWidth',2.2);
hV2 = plot(sig, ann(S.Epi_eu_grid(i_r2)),'--','Color',cEU,'LineWidth',2.2);
set(gca,'YColor','k','XScale','log'); grid on; axis tight;
xlabel('$\sigma_{us}$ (log scale)','Interpreter','latex');
ylabel('Expected inflation, annualized (\%)','Interpreter','latex');
title('(a) $E[\pi''|s]$ by regime','Interpreter','latex');
legend([hU1 hU2 hV1 hV2 hE1 hE2], ...
       {'US, normal','US, volatile','EU, normal','EU, volatile', ...
        'Ergodic, normal','Ergodic, volatile'}, ...
        'Location','best','Interpreter','latex','FontSize',8);

% ---------- (b) Jensen wedge ----------
subplot(1,2,2);
wUS = bpsm(S.Epi_us_grid - S.Epi_us_harm);
wEU = bpsm(S.Epi_eu_grid - S.Epi_eu_harm);
semilogx(sig, wUS(i_r1),'-','Color',cUS,'LineWidth',2.2); hold on;
semilogx(sig, wUS(i_r2),'--','Color',cUS,'LineWidth',2.2);
semilogx(sig, wEU(i_r1),'-','Color',cEU,'LineWidth',2.2);
semilogx(sig, wEU(i_r2),'--','Color',cEU,'LineWidth',2.2);
yline(0,'k:');
grid on; axis tight;
xlabel('$\sigma_{us}$ (log scale)','Interpreter','latex');
ylabel('$E[\Pi''] - (E[1/\Pi''])^{-1}$ (bps, monthly)','Interpreter','latex');
title('(b) Jensen wedge (G07)','Interpreter','latex');
legend({'US, normal','US, volatile','EU, normal','EU, volatile'}, ...
        'Location','best','Interpreter','latex','FontSize',8);

outpng = fullfile('data','expected_inflation_cd.png');
outpdf = fullfile('data','expected_inflation_cd.pdf');
exportgraphics(fig, outpng, 'Resolution', 200);
exportgraphics(fig, outpdf);
fprintf('Saved: %s and %s\n', outpng, outpdf);
