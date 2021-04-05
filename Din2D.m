function dXdt = Din2D(t,X,Foguete)
%DIN2D Função da dinâmica a ser integrada
%Entradas
%               X= r
%                     sigma
%                     u
%                     v
%               Foguete= objeto contendo os dados do foguete

% Define o estágio atual (considera que o tempo inicial começa com a queima
% do primeiro estágio) 
%teste dos estagios



%% Estagio atual
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
%% Dinamica
r=X(1); sigma=X(2); u=X(3); v=X(4); %m=X(5);
%função da massa do foguete, excluindo a massa cc
m=massa(t);

i=Foguete.i;
WE=7.292115856e-5; % modulo da velociodade angular da terra
beta=((pi/2)*(1-t/240)); % CONTROLE
aux=(v-r*WE*cos(i));
gama=atan2(u,aux);
alfa=beta-gama;
RE=6378145;%raio médio da terra.
%% Altitude do Foguete  no tempo (t)
alerta=0;
if r>=RE
    h=r-RE; 
else
    h=0; % caso altitude menor do que zero, zero.
    r=RE;
    alerta=1;
    % caso raio menor do que o raio da terra, raio da terra. % O integrador
    % calcula e armazenta o r anterior ao if, portanto o if mitiga o
    % problema mas o raio ainda pode aparecer como menor que o raio da
    % terra.
end

if t>=173  %
    T=0;%
else%
    T=tracao(h,ie,Foguete); %ie corresponde ao estágio atual
end%
g=AcelGrav(r,0);
% VR é o modulo da velocidaade do veiculo com relação à atmosfera 
VR=(u^2+(v-r*WE*cos(i))^2)^0.5; %i é a inclinação do plano da orbita com relação ao equador.
Mach=VR/Vsom(h);
[FA,FN]=FAero(Mach,alfa,h,VR,ie,Foguete);

rp=u;
sigmap=v/r;
up=(T*sin(beta)-FA*sin(beta)+FN*cos(beta))/m+g+v^2/r;
vp=(T*cos(beta)-FA*cos(beta)-FN*sin(beta))/m-(u*v)/r;
%mp=-Foguete.Mp(ie)/Foguete.tQ(ie); 
if alerta==1
    rp=0; sigmap=0; up=0; vp=0;
end
dXdt=[rp
            sigmap
            up
            vp];
end

