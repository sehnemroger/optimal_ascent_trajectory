function CD=CD0(Mach,ie)
% Função para CD0, atualmente é apenas uma tabela. Ref. Pg. 140.
% Entrada: Mach (numero de mach)
%                ie (estágio atual

%       1 col -> Mach ; 2 col -> CD0
if ie == 1
    %     CDOie=[0, .5, .7, .9, 1.1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5;
    %                  1.9, 1.9, 2.04, 3.07, 4.77, 4.55, 3.79, 3.24, 2.69, 2.45, 2.2, 2.2, 2.2]';
    %     if Mach >= 5
    %         CD=2.2;
    %     else
    %         CD=interp1(CDOie(:,1),CDOie(:,2),Mach);
    % USANDO FUNÇÃO ANALITICA PARA CD
    CD=0.21469E1+(-0.129858E0).*((-1).*(1+exp(1).^(0.326071E2+( ...
        -0.637516E1).*Mach)).^(-1)+(1+exp(1).^((-0.427165E1)+0.124986E2.* ...
        Mach)).^(-1))+(-0.334923E1).*((-1).*(1+exp(1).^((-0.343658E1)+ ...
        0.16902E1.*Mach)).^(-1)+(1+exp(1).^((-0.120308E2)+0.130486E2.* ...
        Mach)).^(-1));
    %     end
    %     CD=2.8462;
elseif ie == 2
    %     CDOie=[0, .5, .7, .9, 1.1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5;
    %         .95, .95, .95, .95, .95, .95, .95, .82, .74, .69, .63, .6, .6]';
    %     if Mach >= 5
    %         CD=.6;
    %     else
    %         CD=interp1(CDOie(:,1),CDOie(:,2),Mach,'pchip');
    % USANDO FUNÇÃO ANALITICA PARA CD
    CD=0.946429E0+(-0.335705E0).*((-1).*(1+exp(1).^(0.122643E2+( ...
        -0.608152E1).*Mach)).^(-1)+(1+exp(1).^(0.694463E1+(-0.188501E1).* ...
        Mach)).^(-1))+(-0.506117E0).*((1+exp(1).^((-0.509511E1)+ ...
        0.821881E0.*Mach)).^(-1)+(-1).*(1+exp(1).^((-0.111885E2)+ ...
        0.520104E1.*Mach)).^(-1));
    %     end
    %     CD=.8254;
else
    CD=0;
end
end