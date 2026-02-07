function [chi_p_psi,thetatemp,psitemp]=Chi_p_psi(mu,p,sigma,iota,lambda,eta)
    F=@(omega,p,sigma) (omega<=0).*exp(p./sigma.*omega).*p+...
        (omega>0).*(p+(1-p).*(1-exp(-(1-p)./sigma.*omega))) ; % Probability short
    Smin=@(mu,p,sigma) sigma.*exp(-p./sigma.*mu)            ;
    Spl =@(mu,p,sigma) mu+sigma.*exp(-p./sigma.*mu)         ; 
    theta=@(mu,p,sigma) Smin(mu,p,sigma)./Spl(mu,p,sigma)   ; % Endogenous Theta 
    LFX_nt_0e_eqs_3;
   %  chi_p=@(theta,iota,lambda) iota*(bartheta(lambda,theta) - bartheta(lambda,theta).^(eta).*theta^(1-eta))/(bartheta(lambda,theta)-1); % Chi Plus formula
    psi  =@(theta,lambda) theta.*(1-exp(-lambda));
    thetatemp = theta(mu,p,sigma);
    psitemp = psi(theta(mu,p,sigma),lambda);
    chi_p_psi=chi_p(thetatemp,iota,lambda,eta)./psitemp;
end