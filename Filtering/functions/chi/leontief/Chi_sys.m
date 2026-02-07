function [if_im,theta_out,psi_out,Smin_out,DW_out,FF_out,Q_out]=Chi_sys(mu,p,sigma,iota,lambda,eta)
    F=@(omega,p,sigma) (omega<=0).*exp(p./sigma.*omega).*p+...
        (omega>0).*(p+(1-p).*(1-exp(-(1-p)./sigma.*omega))) ; % Probability short
    Smin=@(mu,p,sigma) sigma.*exp(-p./sigma.*mu)            ;
    Spl =@(mu,p,sigma) mu+sigma.*exp(-p./sigma.*mu)         ; 
    theta=@(mu,p,sigma) Smin(mu,p,sigma)./Spl(mu,p,sigma)   ; % Endogenous Theta 
    LFX_nt_0e_eqs_3;
    % chi_p=@(theta,iota,lambda) iota*(bartheta(lambda,theta) - bartheta(lambda,theta).^(eta).*theta^(1-eta))/(bartheta(lambda,theta)-1); % Chi Plus formula
    psi  =@(theta,lambda) theta.*(1-exp(-lambda));
    theta_out = theta(mu,p,sigma);
    psi_out = psi(theta(mu,p,sigma),lambda);
    if_im=chi_p(theta(mu,p,sigma),iota,lambda,eta)./psi(theta(mu,p,sigma),lambda);
    Smin_out=Smin(mu,p,sigma);
    Spl_out=Spl(mu,p,sigma);
    DW_out=Smin_out-psi_out*Spl_out;
    FF_out=psi_out*Spl_out;
    Q_out=iota*eta-(chi_p(theta_out,iota,lambda,eta)-chi_m(theta_out,iota,lambda,eta))
end