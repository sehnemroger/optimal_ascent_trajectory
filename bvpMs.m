function sol=bvpMs(ode,jac,bc,solinit,options)
%Attempts to solve a boundary value problem for ODEs, with or without
%unkown parameters, using Multiple Shooting method. [1]
%
%   BVPMS is made to be compatible with MATLAB's bvp5c and bvp4c functions,
%   hence, users can use SOL structure as a initial guess in those functions. More on
%   BVPs see BVP4C reference page.
%
%   options is a strucuture with fields
%       intervals -> Number of divisions on the independent
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
    intervals=12;
    nOuts=1000;
    solver=2; % Choose solver
elseif nargin == 5
    % Options is a structure with those
    intervals=options.intervals;
    nOuts=options.nOuts;
    solver=options.solver;
elseif nargin > 5
    error('Too much entries, see help.');
end

%
nStates=size(solinit.y,1); % number of states in ODE function
nParam=length(solinit.parameters); % Number of parameters
t=linspace(0,1,intervals+1); % This can probably be improved
% based on some information of the initial solution to acelerate
% convergence on problems with continuation.

% Interpolate initial on the shooting points and form the s vector
sInit=reshape(interp1(solinit.x,solinit.y',t(1:end-1),'spline')',[],1);

% Adds additional parameters to the end of the vector
sInit(end+1:end+nParam)=solinit.parameters;

% options for those should be sent with the options structure
if solver == 1
    % Solves the system with fsolve
    FsSolve=@(s)Fs(s,intervals,nStates,nParam,t,ode,bc,jac);
    optsFsolve=optimoptions('fsolve','Display','iter','MaxFunctionEvaluations',1e6,'MaxIterations',1e6,...
        'FunctionTolerance',1e-9,'StepTolerance',1e-9);
    [s,fval]=fsolve(FsSolve,sInit,optsFsolve);
elseif solver == 2
    % Solves the system with fmincon
    [s,fval]=fminconShared(sInit,intervals,nStates,nParam,t,ode,bc,jac);
end
sol=solinit;
sol.parameters=s(end-nParam+1:end); % Parameters
% Recalculates ODEs and accumulates answers
ode1=@(t,X)ode(t,X,sol.parameters);
for i=1:intervals % fazer parfor
    [ti,yi]=ode45(ode1,[t(i) t(i+1)],s(nStates*i-(nStates-1):nStates*i));
    if i==1, y=yi; temp=ti; else, y=[y;yi]; temp=[temp; ti]; end
end
xout=linspace(solinit.x(1),solinit.x(end),nOuts);
[temp,index]=unique(temp);
yout=interp1(temp,y(index,:),xout,'spline')'; % Interpolate final results

% Assemble the solution
sol.x=xout; sol.y=yout; sol.solver='bvp5c'; sol.idata=[]; %those make bpinit give a valid sol structure to bvp5c
sol=bvpinit(sol,[sol.x(1) sol.x(end)]);
sol.stats.maxerr=fval^2;

end %bvpMs

function F=Fs(s,intervals,nStates,nParam,t,odebvp,bc,jac)
tic
param=s(end-nParam+1:end);
uFinal=ones(nStates,10);
ode=@(t,X)odebvp(t,X,param);
% Slice variables for parfor
ti=t; ti(end)=[];
tf=t; tf(1)=[];
s(end)=[];
s=reshape(s,nStates,intervals);
jac1=@(t,X,xp) jac(t,X,[],param);
opts=odeset('Jacobian',jac1);
parfor i=1:intervals
    [~,u]=ode15s(ode,[ti(i) tf(i)],s(:,i),opts); % Solves the ith IVP
    uFinal(:,i)=u(end,:)';
end
uFinal=uFinal(:);
s=s(:);
% Calculates the residual of the continuity equality constraints
internCond=uFinal(1:end-nStates)-s(nStates+1:end);
% Calculates the residual of the contour conditions. See bvp4c bc function.
res=bc(s(1:nStates),uFinal(end-nStates+1:end),param);
% Final Residual
F=[res;internCond];
toc
end %Fs

function [x,f,eflag,outpt]=fminconShared(sInit,intervals,nStates,nParam,t,ode,bc,jac)
FsFmincon1=@(s)FsFmincon(s,intervals,nStates,nParam,t,ode,bc,jac);
con=@(s)Fscon(s,intervals,nStates,nParam,t,ode,bc,jac);
% NESTED FUNCTIONS
% Those are nested to ensure minimum evaluations of the time consuming
% Fs function
sLast=[]; %spam over the nested functions
FsVal=[];
    optsFmincon = optimoptions('fmincon','Display','iter','OptimalityTolerance',1e-9,'StepTolerance',...
        1e-9,'Algorithm','active-set','OutputFcn',@fminconStop);
[x,f,eflag,outpt]=fmincon(FsFmincon1,sInit,[],[],[],[],[],[],con,optsFmincon);
sultimo=sLast;
if FsFmincon1(sultimo)<FsFmincon1(x)
    x=sultimo;
    f=FsFmincon1(sultimo);
end
    function [c,ceq]=Fscon(s,intervals,nStates,nParam,t,odebvp,bc,jac)
        if ~isequal(s,sLast)
            FsVal=Fs(s,intervals,nStates,nParam,t,odebvp,bc,jac);
            sLast=s;
        end
        c=[];
        ceq=FsVal;
    end

    function Fmin=FsFmincon(s,intervals,nStates,nParam,t,odebvp,bc,jac)
        if ~isequal(s,sLast)
            FsVal=Fs(s,intervals,nStates,nParam,t,odebvp,bc,jac);
            sLast=s;
        end
        F=FsVal;
        Fmin=F'*F;
    end
    function stop = fminconStop(s,optimvalues,~)
    % Specifies when fmincon stops
    % That should be and option
    if optimvalues.fval < 1e-10
        stop=true;
        sLast=s; % Catch the out state
    else
        stop=false;
    end
    end
end
