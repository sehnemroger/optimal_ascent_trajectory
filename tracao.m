function T=tracao(h,ie,Foguete)
% Fun��o que calcula a for�a de tra��o. 
% Entradas: (altitude(h), estagio atual(ie), Objeto com os dados do foguete(Foguete))
%                   
% Depende apenas do est�gio, a taxa massica do
% propelente (MPp) � considerada constante ao longo de cada est�gio. (pg.
% 139 da refer�ncia).

As=pi*(Foguete.phie(ie))^2/4;
T=Foguete.CV(ie)*MPp(ie,Foguete)-As*Pa(h);
if T<0
    T=0;
end
% if Foguete.Otimo.Tmax~=0
%     T=Foguete.Otimo.Tmax; % For�a de tra��o constante
% end
end