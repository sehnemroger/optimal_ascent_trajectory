function Mp=MPp(ie,Foguete)
% Taxa de variação massica do propelente em cada estágio considerando
% que a mesma é constante no estágio
% Entradas: ie (estágio)
%                  Foguete (objeto com caracteristicas do foguete em cada
%                                  estágio de voo.

Mp=Foguete.Mp(ie)/Foguete.tQ(ie);

end