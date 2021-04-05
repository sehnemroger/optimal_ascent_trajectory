function [rhr,rhrr]=rho_devs(r)
% Fun��o que calcula a derivada da densidade do ar como fun��o da altitude.
%Modelo de Paulo S�rgio. Pg 156.
%
% Entradas: r [m] (altitude)

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

h=r-6378145; %Raio da terra, por certo deveria vir junto da estrutura Foguete, mas o raio da terra costuma ser uma constante....
h=h/1e3; % Converte h para km
RE=6378145;
if h<=65
    rhr=exp(1).^(((-1)+exp(1).^((-1/1000000000000).*(r+(-1).*RE).*( ...
        1000000000.*k3+(r+(-1).*RE).*(1000000.*k4+(r+(-1).*RE).*(1000.*k5+ ...
        k6.*r+(-1).*k6.*RE))))).*k1+(1/1000).*k2.*((-1).*r+RE)).*(( ...
        -1/1000).*k2+(-1/250000000000).*exp(1).^((-1/1000000000000).*(r+( ...
        -1).*RE).*(1000000000.*k3+(r+(-1).*RE).*(1000000.*k4+(r+(-1).*RE) ...
        .*(1000.*k5+k6.*r+(-1).*k6.*RE)))).*k1.*(250000000.*k3+(r+(-1).* ...
        RE).*(500000.*k4+(r+(-1).*RE).*(750.*k5+k6.*r+(-1).*k6.*RE)))).* ...
        rho0;
    rhrr=(1/62500000000000000000000).*exp(1).^(((-1)+exp(1).^(( ...
        -1/1000000000000).*(r+(-1).*RE).*(1000000000.*k3+(r+(-1).*RE).*( ...
        1000000.*k4+(r+(-1).*RE).*(1000.*k5+k6.*r+(-1).*k6.*RE))))).*k1+( ...
        1/1000).*k2.*((-1).*r+RE)+(-1/500000000000).*(r+(-1).*RE).*( ...
        1000000000.*k3+(r+(-1).*RE).*(1000000.*k4+(r+(-1).*RE).*(1000.*k5+ ...
        k6.*r+(-1).*k6.*RE)))).*((-250000000000).*exp(1).^(( ...
        1/1000000000000).*(r+(-1).*RE).*(1000000000.*k3+(r+(-1).*RE).*( ...
        1000000.*k4+(r+(-1).*RE).*(1000.*k5+k6.*r+(-1).*k6.*RE)))).*k1.*( ...
        500000.*k4+3.*(500.*k5+k6.*(r+(-1).*RE)).*(r+(-1).*RE))+( ...
        250000000.*exp(1).^((1/1000000000000).*(r+(-1).*RE).*(1000000000.* ...
        k3+(r+(-1).*RE).*(1000000.*k4+(r+(-1).*RE).*(1000.*k5+k6.*r+(-1).* ...
        k6.*RE)))).*k2+250000000.*k1.*k3+k1.*(500000.*k4+(750.*k5+k6.*(r+( ...
        -1).*RE)).*(r+(-1).*RE)).*(r+(-1).*RE)).^2+exp(1).^(( ...
        1/1000000000000).*(r+(-1).*RE).*(1000000000.*k3+(r+(-1).*RE).*( ...
        1000000.*k4+(r+(-1).*RE).*(1000.*k5+k6.*r+(-1).*k6.*RE)))).*k1.*( ...
        250000000.*k3+(r+(-1).*RE).*(500000.*k4+(r+(-1).*RE).*(750.*k5+ ...
        k6.*r+(-1).*k6.*RE))).^2).*rho0;
else %elseif h>65% && h<=160
    rhr=0; rhrr=0;
    % %     rhr=exp(1).^((-1).*k1+exp(1).^(((-125)+(1/1000).*(r+(-1).*RE)).*((-1) ...
    % %               .*k8+125.*k9+(-1).*k11.*((-125)+(1/1000).*(r+(-1).*RE)).^3+( ...
    % %               -1/1000).*k9.*(r+(-1).*RE)+(-1).*k10.*(125+(1/1000).*((-1).*r+RE)) ...
    % %               .^2)).*k7+(1/1000).*k2.*((-1).*r+RE)).*((-1/1000).*k2+( ...
    % %               -1/250000000000).*exp(1).^((-1).*((-125)+(1/1000).*(r+(-1).*RE)).* ...
    % %               (k8+k9.*((-125)+(1/1000).*(r+(-1).*RE))+k11.*((-125)+(1/1000).*(r+ ...
    % %               (-1).*RE)).^3+k10.*(125+(1/1000).*((-1).*r+RE)).^2)).*k7.*( ...
    % %               500000.*(500.*k8+k9.*((-125000)+r+(-1).*RE))+k11.*((-125000)+r+( ...
    % %               -1).*RE).^3+750.*k10.*(125000+(-1).*r+RE).^2)).*rho0;
    % %    rhrr=exp(1).^((-1).*k1+exp(1).^(((-125)+(1/1000).*(r+(-1).*RE)).*((-1) ...
    % %               .*k8+125.*k9+(-1).*k11.*((-125)+(1/1000).*(r+(-1).*RE)).^3+( ...
    % %               -1/1000).*k9.*(r+(-1).*RE)+(-1).*k10.*(125+(1/1000).*((-1).*r+RE)) ...
    % %               .^2)).*k7+(1/1000).*k2.*((-1).*r+RE)).*(( ...
    % %               1/62500000000000000000000).*exp(1).^(((-125)+(1/1000).*(r+(-1).* ...
    % %               RE)).*((-1).*k8+125.*k9+(-1).*k11.*((-125)+(1/1000).*(r+(-1).*RE)) ...
    % %               .^3+(-1/1000).*k9.*(r+(-1).*RE)+(-1).*k10.*(125+(1/1000).*((-1).* ...
    % %               r+RE)).^2)).*k7.*(500000.*(500.*k8+k9.*((-125000)+r+(-1).*RE))+ ...
    % %               k11.*((-125000)+r+(-1).*RE).^3+750.*k10.*(125000+(-1).*r+RE).^2) ...
    % %               .^2+(-1/250000000000).*exp(1).^(((-125)+(1/1000).*(r+(-1).*RE)).*( ...
    % %               (-1).*k8+125.*k9+(-1).*k11.*((-125)+(1/1000).*(r+(-1).*RE)).^3+( ...
    % %               -1/1000).*k9.*(r+(-1).*RE)+(-1).*k10.*(125+(1/1000).*((-1).*r+RE)) ...
    % %               .^2)).*k7.*(500000.*k9+1500.*k10.*((-125000)+r+(-1).*RE)+3.*k11.*( ...
    % %               125000+(-1).*r+RE).^2)+((1/1000).*k2+(1/250000000000).*exp(1).^((( ...
    % %               -125)+(1/1000).*(r+(-1).*RE)).*((-1).*k8+125.*k9+(-1).*k11.*(( ...
    % %               -125)+(1/1000).*(r+(-1).*RE)).^3+(-1/1000).*k9.*(r+(-1).*RE)+(-1) ...
    % %               .*k10.*(125+(1/1000).*((-1).*r+RE)).^2)).*k7.*(500000.*(500.*k8+ ...
    % %               k9.*((-125000)+r+(-1).*RE))+k11.*((-125000)+r+(-1).*RE).^3+750.* ...
    % %               k10.*(125000+(-1).*r+RE).^2)).^2).*rho0;
end

end