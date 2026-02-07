function chim = Chi_m(mu,p,sigma,iota,lambda,eta)
    F=@(omega,p,sigma) (omega<=0).*exp(p./sigma.*omega).*p+...
        (omega>0).*(p+(1-p).*(1-exp(-(1-p)./sigma.*omega))) ; % two sided exponential
    Smin=@(mu,p,sigma) sigma.*exp(-p./sigma.*mu)            ; % deficit
    Spl =@(mu,p,sigma) mu+sigma.*exp(-p./sigma.*mu)         ; % surplus
    theta=@(mu,p,sigma) Smin(mu,p,sigma)./Spl(mu,p,sigma)   ; % tightness
    
    MassUnd=@(omega,p,sigma) (omega<=0).*(omega-sigma./p).*exp(p./sigma.*omega).*p+...
       (omega>0).*(-(1-p).*(omega+sigma./(1-p)).*exp(-(1-p)./sigma.*omega));

    bartheta=@(lambda,theta) 1./(1 + (theta.^(-1) - 1) * exp(lambda));
    chi_m=@(theta,iota,lambda,eta) iota.*(bartheta(lambda,theta) - bartheta(lambda,theta).^(eta).*theta.^(-eta))./(bartheta(lambda,theta)-1) ;
    chim=chi_m(theta(mu,p,sigma),iota,lambda,eta);
end