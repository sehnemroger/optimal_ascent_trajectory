function Foguete=perguntas(Foguete)
% Função para perguntas

iinit=input('Alterar condições iniciais? (Enter p/ n) ');
if ~isempty(iinit)
    disp('Para alterar digite o numero. Se desejar não alterar, pressione Enter');
    %% Condições iniciais para orbita inicial circular
    ir0=input(['O raio inicial é  de ',num2str(Foguete.Otimo.r0), ' . Alterando para -> ']);
    if ~isempty(ir0)
        Foguete.Otimo.r0=ir0;
    end
    isigma0=input(['O azimute é  de ',num2str(Foguete.Otimo.sigma0*(pi/180)), ' . Alterando para -> ']);
    if ~isempty(isigma0)
        Foguete.Otimo.sigma0=isigma0/(pi/180);
    end
    Foguete.Otimo.v0=sqrt(Foguete.mue/Foguete.Otimo.r0);
end
ifinal=input('Alterar condições finais? (Enter p/ n) ');
if ~isempty(ifinal)
    disp('Para alterar digite o numero. Se desejar não alterar, pressione Enter');
    %% Condições finais para orbita inicial circular
    irf=input(['O raio final é  de ',num2str(Foguete.Otimo.rf), ' . Alterando para -> ']);
    if ~isempty(irf)
        Foguete.Otimo.rf=irf;
    end
    isigmaf=input(['O azimute final é de ',num2str(Foguete.Otimo.sigmaf*(pi/180)), ' . Alterando para -> ']);
    if ~isempty(isigmaf)
        Foguete.Otimo.sigmaf=isigmaf/(pi/180);
    end
    Foguete.Otimo.vf=sqrt(Foguete.mue/Foguete.Otimo.rf);
end
iinit=input('Alterar fatores de otimalidade? (Enter p/ n) ');
if ~isempty(iinit)
    disp('Para alterar digite o numero. Se desejar não alterar, pressione Enter');
%     %% Fatores de otimilidade 
    iRb=input(['O fator Rb é: ',num2str(Foguete.Otimo.Rb), ' . Alterando para -> ']);
    if ~isempty(iRb)
        Foguete.Otimo.Rb=iRb;
    end
    iR=input(['O fator R é: ',num2str(Foguete.Otimo.R), ' . Alterando para -> ']);
    if ~isempty(iR)
        Foguete.Otimo.R=iR;
    end
    iCT=input(['O fator CT é: ',num2str(Foguete.Otimo.CT), ' . Alterando para -> ']);
    if ~isempty(iCT)
        Foguete.Otimo.CT=iCT;
    end
end