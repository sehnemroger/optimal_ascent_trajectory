function Mp=MPp(ie,Foguete)
% Taxa de varia��o massica do propelente em cada est�gio considerando
% que a mesma � constante no est�gio
% Entradas: ie (est�gio)
%                  Foguete (objeto com caracteristicas do foguete em cada
%                                  est�gio de voo.

Mp=Foguete.Mp(ie)/Foguete.tQ(ie);

end