function Foguete=parametros

%% Parametros
Foguete.phie(1)=0.7; Foguete.phie(2)=1.2; Foguete.phie(3)=0.8; Foguete.phie(4)=0.6; % Diametro da tubeira em cada estágio [m]
Foguete.CV(1)=2550; Foguete.CV(2)=2716; Foguete.CV(3)=2696; Foguete.CV(4)=2762; % Velocidade de exaustão dos gases no vacuo [m/s]
Foguete.Mp(1)=28888; Foguete.Mp(2)=7184; Foguete.Mp(3)=4452; Foguete.Mp(4)=835; % Massa inicial do propelente [kg]
Foguete.tQ(1)=63; Foguete.tQ(2)=62; Foguete.tQ(3)=58; Foguete.tQ(4)=72; % Tempo de queima do propelente [s]
% T1=1.1693e+06; T2=3.1471e+05; T3=2.0694e+05;
Foguete.tQT(1)=63; Foguete.tQT(2)=125; Foguete.tQT(3)=183; Foguete.tQT(4)=255;
Foguete.Sref=0.7964; % area de referencia do foguete [m^2]
Foguete.M(4)=67+136+224; % Massa dos equipamentos, da estrutura e do satélite [kg]
Foguete.M(3)=303+916; % Massa dos equipamentos, da estrutura [kg]
Foguete.M(2)=158+1506; % Massa dos equipamentos, da estrutura [kg]
Foguete.M(1)=5922; % Massa da estrutura [kg]
Foguete.mue=3.986012e14; % m^3 s^-2
%% Parametros adicionais de eficiencia do foguete
Foguete.Otimo.Isp=311; Foguete.Otimo.g0=9.81;