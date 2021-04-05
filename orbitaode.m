function dXdtau=orbitaode(tau,X,region,parametros,Foguete)
% global timeCount timesCounted

% tf1=parametros(1); tf2=parametros(2); tf3=parametros(3);

if region == 1
    t=tau*parametros(1);
    dt_dtau=parametros(1);
else
    t=(tau-region+1)*(parametros(region)-parametros(region-1))+parametros(region-1);
    dt_dtau=parametros(region)-parametros(region-1);
end
% tic
U=controle(t,X,region,Foguete);
% toc
% timeCount=timeCount+toc;
% timesCounted=timesCounted+1;
% if timesCounted >=10000
%     timeCount*1000/10000
%     a=0;
% end
%% Dinamica dos estados
dXdt=Din2D_opt(t,X(1:5),U,region,Foguete);

%% Dinamica dos coestados
r=X(1); sigma=X(2); u=X(3); v=X(4); m=X(5); lr=X(6); ls=X(7); lu=X(8); lv=X(9); lm=X(10);

[FA,FN,FAr,FAu,FAv,FNr,FNu,FNv,FNbeta]=aerodin(t,X(1:5),U,region,Foguete);
[g,gr]=AcelGrav(r); beta=U(1);

switch region
    case 1
        T=Foguete.Otimo.T1;
    case 2
        T=Foguete.Otimo.T2;
    case 3
        T=Foguete.Otimo.T3;        
end


dldt=[(-1).*gr.*lu+r.^(-2).*v.*(ls+(-1).*lv.*u+lu.*v)+FAr.*m.^(-1).*( ...
  lv.*cos(beta)+lu.*sin(beta))+(-1).*FNr.*m.^(-1).*(lu.*cos(beta)+( ...
  -1).*lv.*sin(beta)),0,m.^(-1).*r.^(-1).*((-1).*lr.*m.*r+lv.*m.*v+ ...
  FAu.*r.*(lv.*cos(beta)+lu.*sin(beta))+FNu.*((-1).*lu.*r.*cos(beta) ...
  +lv.*r.*sin(beta))),m.^(-1).*r.^(-1).*((-1).*m.*(ls+(-1).*lv.*u+ ...
  2.*lu.*v)+FAv.*r.*(lv.*cos(beta)+lu.*sin(beta))+FNv.*((-1).*lu.* ...
  r.*cos(beta)+lv.*r.*sin(beta))),m.^(-2).*((-1).*FA.*(lv.*cos(beta) ...
  +lu.*sin(beta))+T.*(lv.*cos(beta)+lu.*sin(beta))+FN.*(lu.*cos( ...
  beta)+(-1).*lv.*sin(beta)))]';

% Dinamica completa
dXdtau=[dXdt
                dldt]*dt_dtau;
end