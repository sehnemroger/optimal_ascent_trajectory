function [g,gr,grr]=AcelGrav(r)
%Aceleração da gravidade. (Modelo do Paulo Sérgio)
% Entradas: r [m] (distancia do centro da terra)

mi=3.986012e14;
g=-mi/r^2;
gr=2.*mi.*r.^(-3);
if nargout >2
    grr=(-6).*mi.*r.^(-4);
end
end