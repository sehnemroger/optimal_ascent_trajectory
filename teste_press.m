% confere se a função da pressão pela altitude está certa 
for i=1:length(h)
alt=h(i)
rh(i)=rho(alt)
P(i)=Pa(h(i))
end
plot(P,h)