function U=controle(t,X,region,Foguete)
% Função que calcula o controle beta de atitude do foguete
% global betai

%
opt=2;
if opt==1
    beta_range=[-pi pi];
    func=@(beta) fun(beta,t,X,region,Foguete);
    gradfunc=@(beta) gradfun(beta,t,X,region,Foguete);
    condc=@(beta) cond(beta,t,X,region,Foguete);
    betai=1.2;
    beta=secTCond(func,gradfunc,condc,beta_range,betai);
%     betai=beta;
elseif opt==2
    beta_range=[-pi pi];
    func=@(beta) fun(beta,t,X,region,Foguete);
    gradfunc=@(beta) gradfun(beta,t,X,region,Foguete);
    grad2func=@(beta) grad2fun(beta,t,X,region,Foguete);
    condc=@(beta) cond(beta,t,X,region,Foguete); 
    betai=1.2;
    beta=step2bCond(func,gradfunc,grad2func,condc,beta_range,betai);
%     betai=beta;
elseif opt==3
    beta=funaprox(X);
end

U=beta;
end

function f=funaprox(X)
lu=X(8); lv=X(9);
f=atan2(-lu,-lv);
end