
clc; clear all;
%% SIMULA O FOGUETE 

%verificar extrapolações do CD
% Definição das váriaveis do foguete
Foguete.phie(1)=0.7; Foguete.phie(2)=1.2; Foguete.phie(3)=0.8; Foguete.phie(4)=0.6; % Diametro da tubeira em cada estágio [m]
Foguete.CV(1)=2550; Foguete.CV(2)=2716; Foguete.CV(3)=2696; Foguete.CV(4)=2762; % Velocidade de exaustão dos gases no vacuo [m/s]
Foguete.Mp(1)=28888; Foguete.Mp(2)=7184; Foguete.Mp(3)=4452; Foguete.Mp(4)=835; % Massa inicial do propelente [kg]
Foguete.tQ(1)=63; Foguete.tQ(2)=62; Foguete.tQ(3)=58; Foguete.tQ(4)=72; % Tempo de queima do propelente [s]
Foguete.tQT(1)=63; Foguete.tQT(2)=125; Foguete.tQT(3)=183; Foguete.tQT(4)=255;
Foguete.Sref=0.7964; % area de referencia do foguete [m^2]
Foguete.M(4)=67+136+224; % Massa dos equipamentos, da estrutura e do satélite [kg]  
Foguete.M(3)=303+916; % Massa dos equipamentos, da estrutura [kg]
Foguete.M(2)=158+1506; % Massa dos equipamentos, da estrutura [kg]
Foguete.M(1)=5922; % Massa da estrutura [kg]
RE=6378145; % raio médio da terra.
WE=7.292115856e-5; % modulo da velociodade angular da terra

%%% CUIDAR COM A MASSA DO PROPELENTE, DIFERENTE DA MASSA DO FOGUETE%%
%%% A MASSA DE PROPELENTE É A MASSA DISPONIVEL PARA CADA ESTAGIO E %%%
%%%%%% % % NÃO A SOMA DE TODOS OS ESTAGIOS RESTANTES%%%%%%%%


%% SIMULAÇÃO


% Condições iniciais
Foguete.i=25*pi/180; 
r0=RE;  % Raio médio da terra, ou seja, lançado da superfície.
sigma0=0;
u0=0; 
v0=465; %WE*RE*cos(i)

X0=[r0 sigma0 u0 v0];

% Integra estagio por estágio considerando que as queimas são sucessivas
% (apenas para teste pois irreal) 
ti=0; tt=[];  Xt=[];
% ESTAGIOS
% for j=1:4
%     [t,X]=ode45(@Din2D,[ti ti+Foguete.tQ(j)],X0,[],j,Foguete);
%     ti=t(end); X0=X(end,:);
%     X0(5)=X0(5)-Foguete.M(j);
%     tt=[tt; t]; Xt=[Xt; X];
% end
tsim=700;
odeoption=odeset('RelTol',1e-9,'AbsTol',1e-9);
[tt,Xt]=ode45(@Din2D,[ti tsim],X0,odeoption,Foguete);

%% GRÁFICOS
% Controle
beta=((pi/2)*(1-tt./240));

% reconstituição dos estágios
vie=zeros(size(tt));
cont=1;
for t=tt'
    if t<=Foguete.tQT(4)
        ie=4;
        if t<=Foguete.tQT(3)
            ie=3;
            if t<=Foguete.tQT(2)
                ie=2;
                if t<=Foguete.tQT(1)
                    ie=1;
                end
            end
        end
    elseif t>Foguete.tQT(4)
        ie=4;
    end
    vie(cont)=ie;
    m(cont)=massa(t);
    cont=cont+1;
end
set(gca,'FontSize',40);
figure(1);
subplot(231); plot(tt,Xt(:,1),'LineWidth',3); xlabel('t [s]','fontsize',40); ylabel('r [m]','fontsize',40);
subplot(232); plot(tt,Xt(:,2),'LineWidth',3); xlabel('t [s]','fontsize',40); ylabel('\sigma [rad]','fontsize',40);
subplot(233); plot(tt,Xt(:,3),'LineWidth',3); xlabel('t [s]','fontsize',40); ylabel('u [m/s]','fontsize',40);
subplot(234); plot(tt,Xt(:,4),'LineWidth',3); xlabel('t [s]','fontsize',40); ylabel('v [m/s]','fontsize',40);
subplot(235); plot(tt,m,'LineWidth',3); xlabel('t [s]','fontsize',40); ylabel('m [kg]','fontsize',40);
subplot(236); plot(tt,beta*180/pi,'LineWidth',3); xlabel('t [s]','fontsize',40); ylabel('\beta [rad]','fontsize',40);