function param=paramAtt2(par,valF,passo)
% Função que calcula o próximo valor do parametro "param" com base
% no valor atual, no passo e no valor final.
%
% Entradas: (par,vals,erro,passoMax)
%       par -> valor atual do parametro
%       valF -> valor final do parametro
%       passo -> passo da variavel
% Saídas: param
%       param -> valor do parametro
%

%% Checa se o valor do parametro reduz ou aumenta
if par < valF
    par_aumenta=true;
elseif par > valF
    par_aumenta=false;
elseif par == valF
    warning('Parametro já se encontra no valor final. ');
    param=par;
    return
end

%% Atualiza o parametro
if par_aumenta
    param=par+passo;
    if param >= valF
        param=valF;
    end
else 
    param=par-passo;
    if param <= valF
        param=valF;
    end
end
end %end function