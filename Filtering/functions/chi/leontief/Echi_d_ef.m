function echi_d = Echi_d_ef(mu,p,sigma,iota,lambda,eta)
    Smin=@(mu,p,sigma) sigma.*exp(-p./sigma.*mu)            ;
    Spl =@(mu,p,sigma) mu+sigma.*exp(-p./sigma.*mu)         ;
    theta=@(mu,p,sigma) Smin(mu,p,sigma)./Spl(mu,p,sigma);
    MassUnd=@(omega,p,sigma) (omega<=0).*(omega-sigma./p).*exp(p./sigma.*omega).*p+...
       (omega>0).*(-(1-p).*(omega+sigma./(1-p)).*exp(-(1-p)./sigma.*omega));
  %  chi_m=@(theta,iota,lambda) iota.*((theta+(1-theta).*exp(lambda)).^(1/2)-theta)./((1-theta).*exp(lambda));
  %  chi_p=@(theta,iota,lambda) iota.*(theta.*(theta+(1-theta).*exp(lambda)).^(1/2)-theta)./((1-theta).*exp(lambda));
     LFX_nt_0e_eqs_3;
  %  theta=@(mu,p,sigma) Smin(mu,p,sigma)./Spl(mu,p,sigma);
  %  chi_p=@(theta,iota,lambda) iota.*(theta.*(theta+(1-theta).*exp(lambda)).^(1/2)-theta)./((1-theta).*exp(lambda));
%     chip=chi_p(theta(mu,p,sigma),iota,lambda,eta);
    echi_d=Echi_d(mu,p,sigma,iota,lambda,eta);
end