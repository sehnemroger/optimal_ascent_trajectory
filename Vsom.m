function Vs=Vsom(h)
% Função para a velocidade do som em função da altitude (Modelo do Paulo
% Sérgio, Pg. 158.)
%
% Entradas: h [m] (altitude)
A0=340;
A1=-7.59813873;
A2=4.41739218e-1;
A3=-1.12391220e-2;
A4=1.70917283e-4;
A5=-1.59888891e-6;
A6=6.51961005e-9;

h=h/1e3; % Converte h para km 

Vs=A0+A1*h+A2*h^2+A3*h^3+A4*h^4+A5*h^5+A6*h^6;

end