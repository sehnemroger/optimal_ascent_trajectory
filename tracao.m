function T=tracao(h,ie,Foguete)
% Função que calcula a força de tração. 
% Entradas: (altitude(h), estagio atual(ie), Objeto com os dados do foguete(Foguete))
%                   
% Depende apenas do estágio, a taxa massica do
% propelente (MPp) é considerada constante ao longo de cada estágio. (pg.
% 139 da referência).

As=pi*(Foguete.phie(ie))^2/4;
T=Foguete.CV(ie)*MPp(ie,Foguete)-As*Pa(h);
if T<0
    T=0;
end
% if Foguete.Otimo.Tmax~=0
%     T=Foguete.Otimo.Tmax; % Força de tração constante
% end
end