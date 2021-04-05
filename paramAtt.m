function param=paramAtt(par,vals,erro,passoMax)
% Função que calcula o próximo valor do parametro "param" com base no
% valor anterior à atualização, no valor atual, no erro causado pela atualização e no valor
% final.
%
% Aplica um controle sigmoidal (similar ao CME) onde a referencia é um
% erro desejado e se controla o passo da variável a ser continuada.
%
% Entradas: (par,vals,erro,passoMax)
%       par -> valor atual do parametro
%       vals -> vetor com valores anterior e final do parametro
%               vals(1) -> valor anterior à atualização
%               vals(2) -> valor final do parametro
%       erro -> erro causado ao bvp4c devido à atualização
%       passoMax -> passo máximo em percentual do valor do parametro
%
% Saídas: param
%       param -> valor do parametro
%
% O valor anterior à atualização define o primeiro incremento, caso 
% não haja ainda valor anterior. Um bom chute pode resultar em melhor (ou pior)
% convergencia.
%
% Roger M. Sehnem

%% Variáveis de controle e parametros de ajuste
er_lim=1e-8; % erro limite, ou de referencia;
pto_conv=1; % valor de erro_er em que o mult é max (aproximadamente);
%% Recolhe valores
valAnt=vals(1);
valF=vals(2);
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
%% Se não há valor anterior (NaN)
if isnan(valAnt)
    if par_aumenta
        valAnt=par-.5*passoMax;
    elseif ~par_aumenta
        valAnt=par+.5*passoMax;
    end
end
%% Cálculo do incremento
incr_old=abs(par-valAnt);  %incremento antigo
erro_er=(er_lim-erro)/er_lim; %erro de referencia do erro
%% controle sigmoidal
sum=-passoMax+passoMax*2*sigmf(erro_er,[pto_conv 0]);
incr=incr_old+sum; %soma ao incremento anterior
if incr >= passoMax
    incr=passoMax;
end
%% Atualiza o parametro
if par_aumenta
    param=par+incr;
    if param >= valF
        param=valF;
    end
else 
    param=par-incr;
    if param <= valF
        param=valF;
    end
end
end %end function