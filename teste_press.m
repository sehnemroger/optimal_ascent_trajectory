% confere se a fun��o da press�o pela altitude est� certa 
for i=1:length(h)
alt=h(i)
rh(i)=rho(alt)
P(i)=Pa(h(i))
end
plot(P,h)