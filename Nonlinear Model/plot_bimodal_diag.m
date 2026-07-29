%% plot_bimodal_diag.m - visualise the bimodality and its cause
N=152; m1=8; m2=4.5;
function [b1,Z1,a1] = slice(mu1,rho1,S1,mu2,rho2,S2,p12,p21,N,m1,m2)
    P=[1-p12 p12; p21 1-p21];
    [Z1,Pi1]=tauchen(N/2,mu1,rho1,S1,m1);
    [Z2,Pi2]=tauchen(N/2,mu2,rho2,S2,m2);
    i12=zeros(N/2,1); i21=zeros(N/2,1);
    for ii=1:N/2
        [~,i12(ii)]=min(abs(Z2-Z1(ii)));
        [~,i21(ii)]=min(abs(Z1-Z2(ii)));
    end
    Pf=zeros(N);
    for ii=1:N/2
        Pf(ii,:)    =[P(1,1)*Pi1(ii,:)      P(1,2)*Pi2(i12(ii),:)];
        Pf(N/2+ii,:)=[P(2,1)*Pi1(i21(ii),:) P(2,2)*Pi2(ii,:)];
    end
    [V,D]=eig(Pf'); [~,jj]=min(abs(diag(D)-1));
    v=real(V(:,jj)); if sum(v)<0, v=-v; end; v=max(v,0); v=v/sum(v);
    b1=v(1:N/2); b1=b1/sum(b1);
    [V1,D1]=eig(Pi1'); [~,j1]=min(abs(diag(D1)-1));
    a1=real(V1(:,j1)); if sum(a1)<0, a1=-a1; end; a1=max(a1,0); a1=a1/sum(a1);
end
rd=@(f) readtable(f);
gC=rd('../Filtering/data/MS_sigma_us_params_C.csv');
gB=rd('../Filtering/data/MS_sigma_us_params_B.csv');
g=@(T,nm) T.value(strcmp(T.param,nm));

fig=figure('Position',[80 80 1040 760],'Color','w');
fmt=@(a) set(a,'FontName','Times','FontSize',12,'Box','on');

% --- Panel A: C, joint slice vs standalone AR(1)
mu1=log(g(gC,'sigma_ss_nor')); r1=g(gC,'rho_nor'); S1=g(gC,'Sigma_nor');
mu2=log(g(gC,'sigma_ss_scr')); r2=g(gC,'rho_scr'); S2=g(gC,'Sigma_scr');
[b1,Z1,a1]=slice(mu1,r1,S1,mu2,r2,S2,g(gC,'trans_nor'),g(gC,'trans_scr'),N,m1,m2);
subplot(2,2,1);
plot(exp(Z1),b1/max(b1),'-','LineWidth',2.6,'Color',[0.85 0.33 0.10]); hold on;
plot(exp(Z1),a1/max(a1),'--','LineWidth',2,'Color',[0.35 0.35 0.35]);
xline(exp(mu2),':','scrambling mean','LineWidth',1.4,'FontSize',9);
grid on; xlim([0.1 0.65]); fmt(gca); xlabel('\sigma_{us}'); ylabel('density (norm.)');
legend('joint-chain slice (model)','standalone AR(1)','Location','northeast');
title(sprintf('C: \\rho_{nor}=%.4f  -> TWO modes', r1));

% --- Panel B: B
mu1b=log(g(gB,'sigma_ss_nor')); r1b=g(gB,'rho_nor'); S1b=g(gB,'Sigma_nor');
mu2b=log(g(gB,'sigma_ss_scr')); r2b=g(gB,'rho_scr'); S2b=g(gB,'Sigma_scr');
[b1b,Z1b,a1b]=slice(mu1b,r1b,S1b,mu2b,r2b,S2b,g(gB,'trans_nor'),g(gB,'trans_scr'),N,m1,m2);
subplot(2,2,2);
plot(exp(Z1b),b1b/max(b1b),'-','LineWidth',2.6,'Color',[0 0.45 0.74]); hold on;
plot(exp(Z1b),a1b/max(a1b),'--','LineWidth',2,'Color',[0.35 0.35 0.35]);
xline(exp(mu2b),':','scrambling mean','LineWidth',1.4,'FontSize',9);
grid on; xlim([0.1 0.65]); fmt(gca); xlabel('\sigma_{us}');
legend('joint-chain slice (model)','standalone AR(1)','Location','northeast');
title(sprintf('B: \\rho_{nor}=%.4f  -> one mode', r1b));

% --- Panel C: rho sweep
subplot(2,2,3);
cols=lines(4); k=0;
for rr=[0.97 0.995 0.9967 0.999]
    k=k+1;
    [bb,ZZ]=slice(mu1,rr,S1,mu2,r2,S2,g(gC,'trans_nor'),g(gC,'trans_scr'),N,m1,m2);
    plot(exp(ZZ),bb/max(bb),'-','LineWidth',2,'Color',cols(k,:)); hold on;
end
grid on; xlim([0.1 0.7]); fmt(gca); xlabel('\sigma_{us}'); ylabel('density (norm.)');
legend('\rho=0.970','\rho=0.995','\rho=0.9967 (C)','\rho=0.999','Location','northeast');
title('Modes multiply as \rho_{nor} \rightarrow 1');

% --- Panel D: half-life vs regime duration
subplot(2,2,4);
rv=linspace(0.90,0.9995,300); hl=log(0.5)./log(rv); dur=1/g(gC,'trans_nor');
semilogy(rv,hl,'-','LineWidth',2.4,'Color',[0.85 0.33 0.10]); hold on;
yline(dur,'--','regime duration (106 mo)','LineWidth',1.8,'FontSize',9);
xline(r1,':','C','LineWidth',1.6,'FontSize',10);
xline(r1b,':','B','LineWidth',1.6,'FontSize',10);
grid on; fmt(gca); xlabel('\rho_{normal}'); ylabel('within-regime half-life (months, log)');
title('Bimodality when half-life > regime duration');

exportgraphics(fig,'output/bimodal_diagnosis.png','Resolution',150);
fprintf('[plot] output/bimodal_diagnosis.png\n');
