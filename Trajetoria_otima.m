% Calculo da trajetória ótima

clc; clear all; close all;

Foguete=parametros; % Puxa parametros do Foguete

%% Condições iniciais para orbita inicial circular
Foguete.i=10*pi/180;
Foguete.Otimo.r0=8.9e6;
Foguete.Otimo.sigma0=0;
Foguete.Otimo.u0=0;
Foguete.Otimo.v0=sqrt(Foguete.mue/Foguete.Otimo.r0);
Foguete.Otimo.m0=5*50591;

%% Condiçoes finais para uma orbita circular
Foguete.Otimo.rf=9e6;
Foguete.Otimo.sigmaf=2*pi; % Apenas para o chute inicial, não é requerido que este seja o sigma final.
Foguete.Otimo.uf=0;
Foguete.Otimo.vf=sqrt(Foguete.mue/Foguete.Otimo.rf);  %velocidade tangencial final para orbita circular
Foguete.Otimo.mf=Foguete.Otimo.m0-200;% Apenas para o chute inicial, não é requerido que esta seja a massa final.

%% Controle Ótimo
% Foguete.Otimo.T=50591;
Foguete.Otimo.k1=10; Foguete.Otimo.k2=10; Foguete.Otimo.k3=10;
Foguete.Otimo.T1=3500000; Foguete.Otimo.T2=3500000; Foguete.Otimo.T3=3500000;
Foguete.Otimo.m12=1e3; Foguete.Otimo.m23=1e3;
Foguete.Otimo.s12=1e3; Foguete.Otimo.s23=1e3;
Foguete.Otimo.Atmosfera=true; % Usar forças aerodinamicas? Deixar em sim

%% Inicialização de Parametros
tolNachieved=true; att_vals=false; att_vals1=true; bvpFail=false;
icont=1; isol=1; contsav=1;

%% Solver
bvpsolver=2; int=50;

%% Teste Controle
% global timeCount timesCounted betai
% timeCount=0; timesCounted=0;
% betai=1.2;
%%
while icont==1
    %% carrega problema anterior
    if isol == 1
        cont_salv=input('Continuar da simulação salva?  (1 p/ sim) ');
        if cont_salv == 1
            load('sol');
            load('Foguete_all');
            Foguete=Foguete_all;
            Foguete.Otimo
            Foguete.i;
            solinit=sol;
            isol=0; tolNachieved=false;
            Foguete.Otimo.Atmosfera=true; % Usar forças aerodinamicas?
        end
    end
    if isol ~= 1 && ~tolNachieved % não é a primeira iteração e a tolerancia foi atingida
        if ~att_vals % entra se flag att_vals é falsa, ou seja, se não esta no modo att automático
            Foguete=perguntas(Foguete); % perguntas sobre mudanças manuais de valores
            att_vals=logical(input('Atualizar os valores automaticamente? (1 p/ s) '));
            if isempty(att_vals), att_vals=false; end
            att_vals1=true;
        end
        if att_vals % modo atualizar valores ativado
            if att_vals1 % caso primeira vez, true por padrão
                wParam=whichParam;
                att_vals1=false;
                valF=input(['Valor final da váriavel ',wParam,' (Atualmente em: ',num2str(Foguete.Otimo.(wParam)),') : ']);
                passo=input('Passo de: (default=.1) '); if isempty(passo), passo=.1; end
            end
            Foguete.Otimo.(wParam)=paramAtt2(Foguete.Otimo.(wParam),valF,passo);
            if Foguete.Otimo.(wParam) == valF
                att_vals=false;
            end
        end
    end
    if isol == 1
        inti=input(['Numero de intervalos do vetor tempo (Aperte enter para utilizar ultimo valor (',num2str(int),')) -> ']);
        if ~isempty(inti), int=inti; end
        t1=linspace(0,1,int/3); t2=linspace(1,2,int/3); t3=linspace(2,3,int/3);
        t=[t1 t2 t3];
        Tf1=100; Tf2=200; Tf3=620; dumblagr1=0; dumblagr2=0;
        parametros(1)=Tf1; parametros(2)=Tf2;  parametros(3)=Tf3;
        parametros(4)=dumblagr1; parametros(5)=dumblagr2;
        inicio1=@(t,k) inicio(t,k,Foguete);
        orbitaode1=@(t,X,region,parametros) orbitaode(t,X,region,parametros,Foguete);
        orbitabc1=@(t,X,parametros) orbitabc(t,X,parametros,Foguete);
        jac1=@(t,X,xp,parametros) jacM(t,X,xp,parametros,Foguete);
        options=bvpset('RelTol',1e-9,'AbsTol',1e-9,'NMax',3e3,'Stats','on');
        solinit=bvpinit(t,inicio1,parametros);
        tic
        if bvpsolver ==1
            soli=bvp4c(orbitaode1,orbitabc1,solinit,options);
        elseif bvpsolver == 2
            soli=bvp5c(orbitaode1,orbitabc1,solinit,options);
        elseif bvpsolver == 3
            soli=mbvpMs(orbitaode1,jac1,orbitabc1,solinit);
        end
        toc
    else
        try
            options=bvpset('RelTol',1e-9,'AbsTol',1e-9,'NMax',3000,'Stats','on');
            orbitaode1=@(t,X,region,parametros) orbitaode(t,X,region,parametros,Foguete);
            orbitabc1=@(t,X,parametros) orbitabc(t,X,parametros,Foguete);
            jac1=@(t,X,xp,parametros) jacM(t,X,xp,parametros,Foguete);
            tic
            if bvpsolver ==1
                soli=bvp4c(orbitaode1,orbitabc1,solinit,options);
            elseif bvpsolver == 2
                soli=bvp5c(orbitaode1,orbitabc1,solinit,options);
            elseif bvpsolver == 3
                soli=mbvpMs(orbitaode1,jac1,orbitabc1,solinit);
            end
            toc
        catch
            warning('O problema resultou em erro. Continuando com a ultima opção salva ')
            load('sol');
            soli=sol;
            bvpFail=true;
        end
    end
    Tf1=soli.parameters(1); Tf2=soli.parameters(2); Tf3=soli.parameters(3);
    X=soli.y; tau=soli.x;
    r=X(1,:); sigma=X(2,:); u=X(3,:); v=X(4,:); m=X(5,:); lr=X(6,:); ls=X(7,:); lu=X(8,:); lv=X(9,:); lm=X(10,:);
%     tau(find(diff(tau)==0))=[];
    firsttau1=true; firsttau2=true;
    for i=1:length(tau)
        if tau(i)<=1 && firsttau1
            t(i)=tau(i)*Tf1;
            if tau(i) == 1, firsttau1=false; end
        elseif tau(i)<=2 && firsttau2
            region=2;
            t(i)=(tau(i)-region+1)*(Tf2-Tf1)+Tf1;
            if tau(i) == 2, firsttau2=false; end
        else
            region=3;
            t(i)=(tau(i)-region+1)*(Tf3-Tf2)+Tf2;
        end
    end
    
    %% Figuras
    % maximo de 4 gráficos por figura para caber de forma aceitável
    cleanfigure
    estados_coestados1=figure(1);
    subplot(2,2,1); plot(t,r,'LineWidth',3); xlabel('t [s]'); ylabel('r [m]','Interpreter','latex');
    subplot(2,2,2); plot(t,lr,'LineWidth',3); xlabel('t [s]'); ylabel('$\lambda_r$','Interpreter','latex');
    subplot(2,2,3); plot(t,sigma/(pi/180),'LineWidth',3); xlabel('t [s]'); ylabel('$\sigma [^\circ/s]$','Interpreter','latex');
    subplot(2,2,4); plot(t,ls,'LineWidth',3); xlabel('t [s]'); ylabel('$\lambda_s$','Interpreter','latex');
    cleanfigure
    estados_coestados2=figure(2);
    subplot(2,2,1); plot(t,u,'LineWidth',3); xlabel('t [s]'); ylabel('u [m/s]','Interpreter','latex');
    subplot(2,2,2); plot(t,lu,'LineWidth',3); xlabel('t [s]'); ylabel('$\lambda_u$','Interpreter','latex');
    subplot(2,2,3); plot(t,v,'LineWidth',3); xlabel('t [s]'); ylabel('v [m/s]','Interpreter','latex');
    subplot(2,2,4); plot(t,lv,'LineWidth',3); xlabel('t [s]'); ylabel('$\lambda_v$','Interpreter','latex');
    cleanfigure
    estados_coestados3=figure(3);
    subplot(1,2,1); plot(t,m,'LineWidth',3); xlabel('t [s]'); ylabel('m [kg]','Interpreter','latex');
    subplot(1,2,2); plot(t,lm,'LineWidth',3); xlabel('t [s]'); ylabel('$\lambda_m$','Interpreter','latex');
    cleanfigure
    
    firstT1=true; firstT2=true;
    for i=1:length(t)
        L=1;
        if t(i)<=Tf1 && firstT1
            Xp=orbitaode(t(i)/Tf1,X(:,i),1,soli.parameters,Foguete);
            U=controle(t(i),X(:,i),1,Foguete);
            H(i)=L+X(6:end,i)'*Xp(1:5)/Tf1;
            if t(i)==Tf1, firstT1=false; end
        elseif t(i)<=Tf2 && firstT2
            Xp=orbitaode(1+(t(i)-Tf1)/(Tf2-Tf1),X(:,i),2,soli.parameters,Foguete);
            U=controle(t(i),X(:,i),2,Foguete);
            H(i)=L+X(6:end,i)'*Xp(1:5)/(Tf2-Tf1);
            if t(i)==Tf2, firstT2=false; end
        else
            Xp=orbitaode(2+(t(i)-Tf2)/(Tf3-Tf2),X(:,i),3,soli.parameters,Foguete);
            U=controle(t(i),X(:,i),3,Foguete);
            H(i)=L+X(6:end,i)'*Xp(1:5)/(Tf3-Tf2);
        end
        phi(i)=U(1);
    end
    
    controles=figure(4);
    subplot(121); plot(t,phi*180/pi,'LineWidth',3); xlabel('t [s]'); ylabel('$\beta\ [^\circ]$','Interpreter','latex');
    subplot(1,2,2); plot(t,H,'LineWidth',3); xlabel('t [s]'); ylabel('H []','Interpreter','latex');
    cleanfigure
    
    %     disp(['Tempo final para inserção em órbita: ', num2str(Tf),' s.'])
    
    % Limpa, por garantia
    clear H u X phi
    
    if soli.stats.maxerr <= 1e-9 && ~att_vals
        sol=soli;
        save('sol.mat','sol');
        Foguete_all=Foguete;
        save('Foguete_all.mat','Foguete_all');
        contsav=contsav+1;
        plot_external_fig=logical(input('Gerar gráficos para LaTeX? (1 p/s Enter p/n)'));
        if plot_external_fig==1
            plot_external_fig=true;
        end
        if plot_external_fig
            matlab2tikz('filename','est_resp1_v12p2.tex','figurehandle',estados_coestados1,'width','0.8\linewidth');%,'extraAxisOptions','ytick distance=3000')
            matlab2tikz('filename','est_resp2_v12p2.tex','figurehandle',estados_coestados2,'width','0.8\linewidth');
            matlab2tikz('filename','est_resp3_v12p2.tex','figurehandle',estados_coestados3,'width','0.8\linewidth');
            matlab2tikz('filename','cont_resp_v12p2.tex','figurehandle',controles,'width','0.8\linewidth');
        end
        iicont=input('Continuar? (Enter p/s 0 p/n)= '); if iicont==0, icont=0; end
        tolNachieved=false;
    end
    if soli.stats.maxerr > 1e-9
        tolNachieved=true;
    else
        tolNachieved=false;
    end
    pause(.0001)
    solinit=soli;
    isol=0;
end