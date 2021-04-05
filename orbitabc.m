function res=orbitabc(Xe,Xd,parametros,Foguete)
%global Foguete

tf1=parametros(1); tf2=parametros(2); tf3=parametros(3);
dumblagr1=parametros(4); dumblagr2=parametros(5);

% Estados iniciais 
ri=Xe(1,1); sigmai=Xe(2,1); ui=Xe(3,1); vi=Xe(4,1); mi=Xe(5,1);
% Estados finais 
r=Xd(1,3); sigma=Xd(2,3); u=Xd(3,3); v=Xd(4,3); m=Xd(5,3);
% Coestados finais
lr=Xd(6,3); ls=Xd(7,3); lu=Xd(8,3); lv=Xd(9,3); lm=Xd(10,3);
% Condições iniciais dos estados
r0=Foguete.Otimo.r0; sigma0=Foguete.Otimo.sigma0; u0=Foguete.Otimo.u0; v0=Foguete.Otimo.v0; m0=Foguete.Otimo.m0;
% Condiçoes finais dos estados
rf=Foguete.Otimo.rf; sigmaf=Foguete.Otimo.sigmaf; uf=Foguete.Otimo.uf; vf=Foguete.Otimo.vf; mue=Foguete.mue;

%Hamiltoniana
Xp1d=orbitaode(1,Xd(:,1),1,parametros,Foguete);
Xp2e=orbitaode(1,Xe(:,2),2,parametros,Foguete);
Xp2d=orbitaode(2,Xd(:,2),2,parametros,Foguete);
Xp3e=orbitaode(2,Xe(:,3),3,parametros,Foguete);
Xp3d=orbitaode(3,Xd(:,3),3,parametros,Foguete);
L=1;
H1d=L+Xd(6:end,1)'*Xp1d(1:5)/tf1;
H2e=L+Xe(6:end,2)'*Xp2e(1:5)/(tf2-tf1);
H2d=L+Xd(6:end,2)'*Xp2d(1:5)/(tf2-tf1);
H3e=L+Xe(6:end,3)'*Xp3e(1:5)/(tf3-tf2);
H3d=L+Xd(6:end,3)'*Xp3d(1:5)/(tf3-tf2);

% massa gasta até a troca do primeiro-segundo estágio
m12=Foguete.Otimo.m12;
% segundo-terceiro
m23=Foguete.Otimo.m23;
% Salto na massa primeiro-segundo
s12=Foguete.Otimo.s12;
% segundo-terceiro
s23=Foguete.Otimo.s23;

res = [ri-r0
        sigmai-sigma0
       ui-u0
       vi-v0
       mi-m0
       %% condicoes primeiro-segundo estagio
       Xd(1:4,1)-Xe(1:4,2)  % estados físicos continuos, exceto a massa
       Xd(6:9,1)-Xe(6:9,2) % Continuidade dos coestados, mais no mathematica
       Xe(10,2)-dumblagr1
       Xd(5,1)-(m0-m12)
       Xd(5,1)-s12-Xe(5,2)
       %% condicoes segundo-terceiro estagio
       Xd(1:4,2)-Xe(1:4,3)  % estados físicos continuos, exceto a massa
       Xd(6:9,2)-Xe(6:9,3) % Continuidade dos coestados, mais no mathematica
       Xe(10,3)-dumblagr2
       Xd(5,2)-(m0-m12-s12-m23)
       Xd(5,2)-s23-Xe(5,3)
       r-rf
       ls
       u
       v-sqrt(mue/r)
       lm
       H1d-H2e
       H2d-H3e
       H3d];
end