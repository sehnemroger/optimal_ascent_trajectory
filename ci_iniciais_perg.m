function Foguete=ci_iniciais_perg(Foguete)
% Fun��o para perguntas

fields=fieldnames(Foguete.Otimo);
otimo=struct2cell(Foguete.Otimo);
for i=1:4
    in=input([fields{i},' � ',num2str(otimo{i}), ' . Alterando para -> ']);
    if ~isempty(in)
        otimo{i}=in;
    end
end

Foguete.Otimo=cell2struct(otimo,fields);