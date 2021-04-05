function g2f=grad2fun(beta,t,X,region,Foguete)

r=X(1); sigma=X(2); u=X(3); v=X(4); m=X(5); lr=X(6); ls=X(7); lu=X(8); lv=X(9); lm=X(10);

U(1)=beta;
[FA,FN,FAr,FAu,FAv,FNr,FNu,FNv,FNbeta]=aerodin(t,X,U,region,Foguete);

switch region
    case 1
        T=Foguete.Otimo.T1;
    case 2
        T=Foguete.Otimo.T2;
    case 3
        T=Foguete.Otimo.T3;        
end

g2f=m.^(-1).*(FN.*(lv.*cos(beta)+lu.*sin(beta))+FA.*(lu.*cos(beta)+( ...
  -1).*lv.*sin(beta))+(-1).*(3.*FNbeta+T).*(lu.*cos(beta)+(-1).*lv.* ...
  sin(beta)));
end