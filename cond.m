function c=cond(beta,t,X,region,Foguete)

%% restriçoes de inequalidade c<=0
autoValores=eigHuu(beta,t,X,region,Foguete);
c=[-autoValores,beta-pi,-beta-pi];
end