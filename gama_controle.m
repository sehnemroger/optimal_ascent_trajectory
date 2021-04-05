function [gamac,ggamac,g2gamac]=gama_controle(beta,t,X,Foguete)

r=X(1); sigma=X(2); u=X(3); v=X(4); m=X(5); lr=X(6); ls=X(7); lu=X(8); lv=X(9); lm=X(10);

g0=Foguete.Otimo.g0; Isp=Foguete.Otimo.Isp; R=Foguete.Otimo.R; Tmax=Foguete.Otimo.Tmax;
gamac=(1/2).*g0.^(-1).*Isp.^(-1).*m.^(-1).*R.^(-1).*Tmax.*(lm.*m+(-1).* ...
  g0.*Isp.*(lv.*cos(beta)+lu.*sin(beta)));

if nargout > 1
    ggamac=(-1/2).*m.^(-1).*R.^(-1).*Tmax.*(lu.*cos(beta)+(-1).*lv.*sin(beta));
    if nargout > 2
        g2gamac=(-1/2).*m.^(-1).*R.^(-1).*Tmax.*((-1).*lv.*cos(beta)+(-1).*lu.* ...
                          sin(beta));
    end
end
end