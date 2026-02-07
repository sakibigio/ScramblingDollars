function echi_m = Echi_m(mu,p,sigma,iota,lambda,eta)
    F=@(omega,p,sigma) (omega<=0).*exp(p./sigma.*omega).*p+...
        (omega>0).*(p+(1-p).*(1-exp(-(1-p)./sigma.*omega)));
    Smin=@(mu,p,sigma) sigma.*exp(-p./sigma.*mu)            ; % deficit
    Spl =@(mu,p,sigma) mu+sigma.*exp(-p./sigma.*mu)         ; % surplus
    MassUnd=@(omega,p,sigma) (omega<=0).*(omega-sigma./p).*exp(p./sigma.*omega).*p+...
       (omega>0).*(-(1-p).*(omega+sigma./(1-p)).*exp(-(1-p)./sigma.*omega));
    theta=@(mu,p,sigma) Smin(mu,p,sigma)./Spl(mu,p,sigma)   ; % tightness
    bartheta=@(lambda,theta) 1./(1 + (theta^(-1) - 1) * exp(lambda));
    chi_m=@(theta,iota,lambda,eta) iota.*(bartheta(lambda,theta) - bartheta(lambda,theta).^(eta).*theta^(-eta))/(bartheta(lambda,theta)-1) ;
    chi_p=@(theta,iota,lambda,eta) iota.*(bartheta(lambda,theta) - bartheta(lambda,theta).^(eta).*theta^(1-eta))/(bartheta(lambda,theta)-1);

    % Marginal Liquidity at Interior
    Echi_m=@(mu,p,sigma,iota,lambda,eta) (1-F(-mu,p,sigma)).*chi_p(theta(mu,p,sigma),iota,lambda,eta)...
        +F(-mu,p,sigma).*chi_m(theta(mu,p,sigma),iota,lambda,eta);
    
    echi_m=Echi_m(mu,p,sigma,iota,lambda,eta);
end
