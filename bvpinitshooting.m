function solinit=bvpinitshooting(ode,odebc,initbc,int,Foguete)
% Função para encontar uma solução inicial atravéz
% do método dos multiplos tiros, parametros desconhecidos
%  apenas no tempo
%
% initbc contém os chutes dos valores de contorno iniciais em y(0) e
% do tempo inicial
% initbc=[y[0]' ti]

func=@(bc) fun(bc,ode,odebc);
ipattern=input('Usar Patternsearch (=1), ga (=2) ou surrogate (=3) --> ');
if ipattern == 1
    options = optimoptions('patternsearch','Display','iter','PlotFcn',@psplotbestf,'UseParallel',true,'UseCompletePoll'...
        ,true,'MeshContractionFactor',.5,'Cache','on','CacheSize',1e4,'MaxFunctionEvaluations',1e7,'MaxIterations',1e7,'SearchFcn',...
     'MADSPositiveBasisNp1');
    initbc=patternsearch(func,initbc,[],[],[],[],[],[],[],options);
elseif ipattern == 2
    options = optimoptions('ga','Display','iter','PlotFcn',{@gaplotdistance,...
        @gaplotscores,@gaplotbestf,@gaplotbestindiv},...
        'PopulationSize',100,'UseParallel',true);
    initbc=ga(func,length(initbc),[],[],[],[],[],[],[],options);
elseif ipattern == 3
    ub=[Foguete.Otimo.rf .2*pi 1e5 1e10 1e20 1e20  1e20  1e20  100000];
    lb=[Foguete.Otimo.r0 0.0000001 0.000001 0.000001 -1e20  -1e20  -1e20  -1e20  .1];
    options = optimoptions('surrogateopt','PlotFcn','surrogateoptplot',...
    'MaxFunctionEvaluations',1e8,'UseParallel',true,'ObjectiveLimit',1e-5,'Display','iter','InitialPoints',initbc);
    initbc=surrogateopt(func,lb,ub,options);
end
odec=@(t,X) ode(t,X,initbc(end));
sol=ode45(odec,[0 initbc(end)],initbc(1:end-1));
param=initbc(end);
t=linspace(0,1*param,int);
X=deval(sol,t);
solinit.parameters=param;
solinit.x=t; solinit.y=X;
%% plots
r=X(1,:); sigma=X(2,:); u=X(3,:); v=X(4,:); lr=X(5,:); ls=X(6,:); lu=X(7,:); lv=X(8,:);
figure(1)
subplot(241); plot(t,r,'LineWidth',3); xlabel('t [s]','fontsize',40); ylabel('r [m]','Interpreter','latex','fontsize',40);
subplot(242); plot(t,sigma/(pi/180),'LineWidth',3); xlabel('t [s]','fontsize',40); ylabel('$\sigma [^\circ/s]$','Interpreter','latex','fontsize',40);
subplot(243); plot(t,u,'LineWidth',3); xlabel('t [s]','fontsize',40); ylabel('u [m/s]','Interpreter','latex','fontsize',40);
subplot(244); plot(t,v,'LineWidth',3); xlabel('t [s]','fontsize',40); ylabel('v [m/s]','Interpreter','latex','fontsize',40);
subplot(245); plot(t,lr,'LineWidth',3); xlabel('t [s]','fontsize',40); ylabel('$\lambda_r$','Interpreter','latex','fontsize',40);
subplot(246); plot(t,ls,'LineWidth',3); xlabel('t [s]','fontsize',40); ylabel('$\lambda_s$','Interpreter','latex','fontsize',40);
subplot(247); plot(t,lu,'LineWidth',3); xlabel('t [s]','fontsize',40); ylabel('$\lambda_u$','Interpreter','latex','fontsize',40);
subplot(248); plot(t,lv,'LineWidth',3); xlabel('t [s]','fontsize',40); ylabel('$\lambda_v$','Interpreter','latex','fontsize',40);
end

function res=fun(bc,ode,odebc)

odec=@(t,X) ode(t,X,bc(end));
[~,X]=ode45(odec,[0 bc(end)],bc(1:end-1));
Xb=X(end,:); Xa=X(1,:);
res=sum(abs(odebc(Xa,Xb,bc(end))));
end