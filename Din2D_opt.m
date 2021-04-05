function dXdt = Din2D_opt(t,X,U,region,Foguete)
%DIN2D Função da dinâmica a ser integrada
%Entradas
%               X= r
%                     sigma
%                     u
%                     v
%               Foguete= objeto contendo os dados do foguete
%               u= beta -> controle (atitude do foguete)
%                     gama -> controle de tração
% Define o estágio atual (considera que o tempo inicial começa com a queima
% do primeiro estágio)

%% Dinamica
r=X(1); sigma=X(2); u=X(3); v=X(4); m=X(5);
beta=U(1);
[FA,FN]=aerodin(t,X,U,region,Foguete);

g=AcelGrav(r);

switch region
    case 1
        ki=Foguete.Otimo.k1;
        T=Foguete.Otimo.T1;
    case 2
        ki=Foguete.Otimo.k2;
        T=Foguete.Otimo.T2;
    case 3
        ki=Foguete.Otimo.k3;
        T=Foguete.Otimo.T3;
end
% rp=u;
% sigmap=v/r;
% up=(T*sin(beta)-FA*sin(beta)+FN*cos(beta))/m+g+v^2/r;
% vp=(T*cos(beta)-FA*cos(beta)-FN*sin(beta))/m-(u*v)/r;
% mp=-ki;

dXdt=[u,r.^(-1).*v,g+r.^(-1).*v.^2+m.^(-1).*(FN.*cos(beta)+(-1).*FA.* ...
  sin(beta)+T.*sin(beta)),(-1).*r.^(-1).*u.*v+m.^(-1).*((-1).*FA.* ...
  cos(beta)+T.*cos(beta)+(-1).*FN.*sin(beta)),(-1).*ki]';
end

