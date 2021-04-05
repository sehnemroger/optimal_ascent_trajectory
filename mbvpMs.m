function sol=mbvpMs(ode,jac,bc,solinit,options)
%Attempts to solve a boundary value problem for ODEs, with or without
%unkown parameters, using Multiple Shooting method. [1]
%
%   BVPMS is made to be compatible with MATLAB's bvp5c and bvp4c functions,
%   hence, users can use SOL structure as a initial guess in those functions. More on
%   BVPs see BVP4C reference page.
%
%   options is a structure with fields
%       intervals -> Number of divisions per region on the independent
%       variable to solve the problem.
%       nOuts -> Number of desired outputs per state.
%       solver -> Solver for the resulting nonlinear system of
%       equations.
%           1 -> Uses fsolve;
%           2 -> Uses fmincon;
%           
%           Note: fmincon tends to use less functions evaluations.
%
% References:
%       [1]  Keskin AÜ. Multiple Shooting Method. In: Boundary Value
%              Problems for Engineers

% Róger M. Sehnem, 2020.


% Just work when all the options structure is complete, improove that.
if nargin < 4
    error('Few entries, see help.');
elseif nargin < 5
    % No options, using default ones.
    iintervals=2;
    nOuts=3000;
    solver=2; % Choose solver
elseif nargin == 5
    % Options is a structure with those
    iintervals=options.iintervals;
    nOuts=options.nOuts;
    solver=options.solver;
elseif nargin > 5
    error('Too much entries, see help.');
end

%
nStates=size(solinit.y,1); % number of states in ODE function
nParam=length(solinit.parameters); % Number of parameters
%%
intbci=find(diff(solinit.x)==0); % localiza pontos de interface
nreg=length(intbci)+1; % number of regions
t=zeros(nreg,iintervals+1);
for i=1:nreg
    if i==1
        t(i,:)=linspace(solinit.x(1),solinit.x(intbci(i)),iintervals+1); % rounds to the closest lower value
    elseif i~=nreg
        t(i,:)=linspace(solinit.x(intbci(i-1)),solinit.x(intbci(i)),iintervals+1);
    else
        t(i,:)=linspace(solinit.x(intbci(i-1)),solinit.x(end),iintervals+1);
    end
end

t=reshape(t',[],1); % transform into a column vector
% solinit.x(intbci)=[];  solinit.y(:,intbci)=[];
intbc=find(diff(t)==0); % localiza pontos de interface no novo tempo
tintbc=t(intbc); %pega os valores das discontinuidades
t(intbc)=[]; %limpa o tempo repetido
for i=1:length(tintbc)
    intbc(i)=find(t==tintbc(i));
end
intervals=length(t)-1;
%%
% t=linspace(0,1,intervals+1);
% Iterpolate initial guess by region
sInit=zeros(nStates,intervals);
aux=0;
for i=1:nreg
    if i==1
        sInit(:,1:iintervals)=interp1(solinit.x(1:intbci(i)),solinit.y(:,1:intbci(i))',t(1:iintervals),'pchip')';
    elseif i==nreg
        sInit(:,aux+1:aux+iintervals)=interp1(solinit.x(intbci(i-1)+1:end),solinit.y(:,intbci(i-1)+1:end)',t(aux+1:aux+iintervals),'pchip')';
    else
        sInit(:,aux+1:aux+iintervals)=interp1(solinit.x(intbci(i-1)+1:intbci(i)),solinit.y(:,intbci(i-1)+1:intbci(i))',t(aux+1:aux+iintervals),'pchip')';
    end
    aux=aux+iintervals;
end

sInit=reshape(sInit,[],1);

% Interpolate initial on the shooting points and form the s vector
% sInit=reshape(interp1(solinit.x,solinit.y',t(1:end-1))',[],1);


% Adds additional parameters to the end of the vector
sInit(end+1:end+nParam)=solinit.parameters;

% options for those should be sent with the options structure
if solver == 1
    % Solves the system with fsolve
    FsSolve=@(s)Fs(s,intervals,nreg,tintbc,intbc,nStates,nParam,t,ode,bc,jac);
    optsFsolve=optimoptions('fsolve','Display','iter','MaxFunctionEvaluations',1e6,'MaxIterations',1e6,...
        'FunctionTolerance',1e-9,'StepTolerance',1e-9,'Algorithm','trust-region-dogleg');
    [s,fval]=fsolve(FsSolve,sInit,optsFsolve);
elseif solver == 2
    % Solves the system with fmincon
    [s,fval]=fminconShared(sInit,intervals,nreg,tintbc,intbc,nStates,nParam,t,ode,bc,jac);
end
sol=solinit;
sol.parameters=s(end-nParam+1:end); % Parameters
% Recalculates ODEs and accumulates answers
for i=1:intervals
    region=1;
    for j=1:nreg-1
        if t(i)>=tintbc(j)
            region=region+1;
        end
    end
    oder{i}=@(t,X)ode(t,X,region,sol.parameters);
end
jac1=@(t,X,xp) jac(t,X,[],sol.parameters);
opts=odeset('Jacobian',jac1,'Refine',1);
for i=1:intervals % fazer parfor
    [ti,yi]=ode15s(oder{i},[t(i),t(i+1)],s(nStates*i-(nStates-1):nStates*i),opts);
    if i==1, y=yi; temp=ti; else, y=[y;yi]; temp=[temp; ti]; end
end
for i=1:nreg
    if i==1
        xout(i,:)=linspace(solinit.x(1),tintbc(1),nOuts); 
    elseif i~=nreg
        xout(i,:)=linspace(tintbc(i-1),tintbc(i),nOuts);
    else
        xout(i,:)=linspace(tintbc(end),solinit.x(end),nOuts);
    end
end
xout=reshape(xout',[],1); % transform into a column vector
% interp separadamente
 % Interpolate final results
intbc=zeros(2,length(tintbc));
for i=1:length(tintbc)
    intbc(:,i)=find(temp==tintbc(i)); % localiza pontos de interface no novo temp
end
 intbc=reshape(intbc,[],1);
% intbc=find(temp==tintbc); % localiza pontos de interface no novo temp
%%
yout=zeros(nStates,length(xout));
aux=0; aux2=0;
for i=1:nreg
    if i==1
        [t,index]=unique(temp(1:intbc(1)));
        yu=y(index,:)';
        yout(:,1:nOuts)=interp1(t,yu',xout(1:nOuts),'pchip')';
    elseif i==nreg
        [t,index]=unique(temp(intbc(aux2):end));
        index=index+intbc(aux2)-1;
        yu=y(index,:)';
        yout(:,aux+1:aux+nOuts)=interp1(t,yu',xout(aux+1:aux+nOuts),'pchip')';
    else
        [t,index]=unique(temp(intbc(aux2):intbc(aux2+1)));
        index=index+intbc(aux2)-1;
        yu=y(index,:)';
        yout(:,aux+1:aux+nOuts)=interp1(t,yu',xout(aux+1:aux+nOuts),'pchip')';
    end
    aux=aux+nOuts;
    aux2=aux2+2;
end

% Assemble the solution
sol.x=xout'; sol.y=yout; sol.solver='bvp5c'; sol.idata=[]; %those make bpinit give a valid sol structure to bvp5c
sol=bvpinit(sol,[sol.x(1) sol.x(end)]);
sol.stats.maxerr=fval^2;

end %bvpMs

function F=Fs(s,intervals,nreg,tintbc,intbc,nStates,nParam,t,odebvp,bc,jac)
% tic
param=s(end-nParam+1:end);
uFinal=ones(nStates,intervals);
% ode=@(t,X)odebvp(t,X,param);
% Slice variables for parfor
ti=t; ti(end)=[];
tf=t; tf(1)=[];
for i=1:intervals
    region=1;
    for j=1:nreg-1
        if ti(i)>=tintbc(j)
            region=region+1;
        end
    end
    ode{i}=@(t,X)odebvp(t,X,region,param);
end
s(end-nParam+1:end)=[]; % apaga os parametros
s=reshape(s,nStates,intervals);
jac1=@(t,X,xp) jac(t,X,[],param);
opts=odeset('Jacobian',jac1);
parfor i=1:intervals %parfor
    [~,u]=ode15s(ode{i},[ti(i) tf(i)],s(:,i),opts); % Solves the ith IVP
    uFinal(:,i)=u(end,:)';
end
interf1=1; interf2=1;
for i=1:nreg
    if i==1
        Xe(:,i)=s(:,1);
        rmcols(i)=1;
    else
        Xe(:,i)=s(:,intbc(interf1));
        rmcols(i)=intbc(interf1);
        interf1=interf1+1;
    end
    
    if i==nreg
        Xd(:,i)=uFinal(:,end);
        rmcoluf(i)=size(uFinal,2);
%         uFinal(:,end)=[];
    else
        Xd(:,i)=uFinal(:,intbc(interf2)-1);
        rmcoluf(i)=intbc(interf2)-1;
%         uFinal(:,intbc(interf2)-1)=[];
        interf2=interf2+1;
    end
end
uFinal(:,rmcoluf)=[];
s(:,rmcols)=[];
uFinal=uFinal(:);
s=s(:);
% Calculates the residual of the continuity equality constraints
% internCond=uFinal(1:end-nStates)-s(nStates+1:end);
internCond=uFinal-s; % Lindamente simples
% Calculates the residual of the contour conditions. See bvp4c bc function.
res=bc(Xe,Xd,param);
% Final Residual
F=[res;internCond];
% toc
end %Fs

function [x,f,eflag,outpt]=fminconShared(sInit,intervals,nreg,tintbc,intbc,nStates,nParam,t,ode,bc,jac)
FsFmincon1=@(s)FsFmincon(s,intervals,nreg,tintbc,intbc,nStates,nParam,t,ode,bc,jac);
con=@(s)Fscon(s,intervals,nreg,tintbc,intbc,nStates,nParam,t,ode,bc,jac);
% NESTED FUNCTIONS
% Those are nested to ensure minimum evaluations of the time consuming
% Fs function
sLast=[]; %spam over the nested functions
FsVal=[];
    optsFmincon = optimoptions('fmincon','Display','iter','OptimalityTolerance',1e-12,'StepTolerance',...
        1e-12,'Algorithm','active-set','OutputFcn',@fminconStop);
[x,f,eflag,outpt]=fmincon(FsFmincon1,sInit,[],[],[],[],[],[],con,optsFmincon);
sultimo=sLast;
if FsFmincon1(sultimo)<FsFmincon1(x)
    x=sultimo;
    f=FsFmincon1(sultimo);
end
    function [c,ceq]=Fscon(s,intervals,nreg,tintbc,intbc,nStates,nParam,t,odebvp,bc,jac)
        if ~isequal(s,sLast)
            FsVal=Fs(s,intervals,nreg,tintbc,intbc,nStates,nParam,t,odebvp,bc,jac);
            sLast=s;
        end
        c=[];
        ceq=FsVal;
    end

    function Fmin=FsFmincon(s,intervals,nreg,tintbc,intbc,nStates,nParam,t,odebvp,bc,jac)
        if ~isequal(s,sLast)
            FsVal=Fs(s,intervals,nreg,tintbc,intbc,nStates,nParam,t,odebvp,bc,jac);
            sLast=s;
        end
        F=FsVal;
        Fmin=F'*F;
    end
    function stop = fminconStop(s,optimvalues,~)
    % Specifies when fmincon stops
    % That should be and option
    if optimvalues.fval < 1e-12
        stop=true;
        sLast=s; % Catch the quiting state
    else
        stop=false;
    end
    end
end
