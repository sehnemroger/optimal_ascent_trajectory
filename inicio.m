function f = inicio(tau,k,Foguete)

% Condições iniciais dos estados
r0=Foguete.Otimo.r0; sigma0=Foguete.Otimo.sigma0;
u0=Foguete.Otimo.u0; v0=Foguete.Otimo.v0; 
m0=Foguete.Otimo.m0;
% Condiçoes finais dos estados
rf=Foguete.Otimo.rf; sigmaf=Foguete.Otimo.sigmaf; 
uf=Foguete.Otimo.uf; vf=Foguete.Otimo.vf; 
mf=Foguete.Otimo.mf;
f=[r0*(1-tau)+tau*rf
       sigma0*(1-tau)+tau*sigmaf
       u0*(1-tau)+tau*uf
       v0*(1-tau)+tau*vf
       m0*(1-tau)+tau*mf
       1
       0
       1
       1
       1];
end