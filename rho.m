function rh=rho(h)
% Função que calcula a densidade do ar como função da altitude (Paulo
% Sérgio. Pg 156
%
% Entradas: h [m] (altitude)

rho0=1.225;
k1=1.02280550;
k2=1.21226930e-1;
k3=-3.48643241e-2;
k4=3.50991865e-3;
k5=-8.33000533e-5;
k6=1.15219733e-6;
k7=2.204504;
k8=3.12550644e-3;
k9=5.82501032e-4;
k10=7.44266792e-6;
k11=3.88398609e-7;

h=h/1e3; % Converte h para km 

if h<=65
    phi=k1*exp(-(k3*h+k4*h^2+k5*h^3+k6*h^4));
    rh=rho0*exp(-(k1+k2*h-phi));
else %elseif h>65 && h<=160
    rh=0;
%     csi=h-125;
%     psi=k7*exp(-(k8*csi+k9*csi^2+k10*csi^3+k11*csi^4));
%     rh=rho0*exp(-(k1+k2*h-psi));
end
end