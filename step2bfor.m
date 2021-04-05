function s=step2b(fun,grad,grad2,s)
% Finds s when F(s)=0 based on the Two Step Iterative Method
% for Finding Root of a Nonlinear Equation². Use analytical
% gradients for improved convergence rate.
%
% Inputs: fun           =      F(s) {function handle}
%              grad          =     dF(s)/ds {function gandle}
%              grad2        =     d²F(s)/d²s {function gandle}
%              s                =      initial gues of F(s)=0
% Outputs: s
%
% When using analytical derivatives the algorithm uses a maximum value for
% them, this improve convergence. Source: by experience.
%
% Can be used without providing analytical derivatives. In that case,
% the algotithm wont use any possible derivatives. Meaning that
% s=step2b(fun,grad,[],s) is the same as s=step2b(fun,[],[],s).
%
% When using the algorithm without the analytical derivatives, there is
% another "option", cMax, that controls how long should the algorithm use
% low order (less expensive) derivatives.
%
% Those options can be easily implemented, but, by now, they're hard
% coded options.

% Method by Rajesh C. Shah, Rajiv B. Shah in "Two Step Iterative Method
% for Finding Root of a Nonlinear Equation".
%
% Implemented by Róger M. Sehnem, 2020.

%% Decides if it's using analytical derivatives or not
if isempty(grad) || isempty(grad2)
    analytic=false;
else
    analytic=true;
end
%% Main loop
counter=0; % initialize the counter
max_error=1e-9; % should be an option

if analytic
    for i=1:10
        s = twoStepAnalytic(fun,grad,grad2,s);
        counter=counter+1;
    end
else
    while abs(fun(s)) >= max_error
    s = twoStepFiniteDiff(fun,s,counter);
    counter=counter+1;
    end
end
% counter
end

function s=twoStepAnalytic(f,grad,grad2,s)
%% Calculate function values and derivatives
feval=f(s);
gradf=grad(s);
gradf2=grad2(s);
%% set maximum absolute value for derivatives
gradmax=150; % -> 2
gradmin=-gradmax;
gradf=min(max(gradf,gradmin),gradmax);
gradf2=min(max(gradf2,gradmin),gradmax);
%
m=2; % Works with that but can be larger and should be an option.
%
%% Calculate next step (1st step)
rootPart=sqrt(1/m^2*feval^2*gradf^2-(1-1/m^2)*feval^3*gradf2);
if ~isreal(rootPart), rootPart=0; end
%                            - -> 1
yn=s-(feval*gradf+rootPart)/(feval*gradf2+gradf^2);
% Calculate next step (2nd step)
fevalyn=f(yn);
gradfyn=grad(yn);
% gradfyn=min(max(gradfyn,gradmin),gradmax);
s=yn-fevalyn/gradfyn;
end

function phi=twoStepFiniteDiff(f,phi,counter)
% Function variables are have not  the same name because i'm too lazy
delt=.5e-8;
feval=f(phi);
cMax=50; % -> 3
%
if counter <= cMax
    % Low order estimate
    fma1=f(phi+1*delt); fma2=f(phi+2*delt);
    gradf=(fma1-feval)/(delt);
    gradf2=(feval-2*fma1+fma2)/(delt^2);
else
    % High order estimate
    fme2=f(phi-2*delt); fme1=f(phi-1*delt); fma1=f(phi+1*delt); fma2=f(phi+2*delt);
    gradf=(fme2-8*fme1+8*fma1-fma2)/(12*delt);
    gradf2=-(fme2-16*fme1+30*feval-16*fma1+fma2)/(12*delt^2);
end
%
m=2;
%
rootPart=sqrt(1/m^2*feval^2*gradf^2-(1-1/m^2)*feval^3*gradf2);
if ~isreal(rootPart), rootPart=0; end
%                               -
yn=phi-(feval*gradf+rootPart)/(feval*gradf2+gradf^2);
fevalyn=f(yn);
if counter <= cMax
    % Low order estimate
    fma1yn=f(yn+1*delt);
    gradfyn=(fma1yn-fevalyn)/(delt);
else
    % High order estimate
    fme2yn=f(yn-2*delt); fme1yn=f(yn-1*delt); fma1yn=f(yn+1*delt); fma2yn=f(yn+2*delt);
    gradfyn=(fme2yn-8*fme1yn+8*fma1yn-fma2yn)/(12*delt);
end
phi=yn-fevalyn/gradfyn;
end

%% Comments
% 1 -> That part can be + or -,
% i put it as only a + but, ideally, i think that it should
% be changed automactly if it get stucked (it
% happends, sometimes, that the algorithm get
% stuck between 2, 3 or more points).
%
% 2 -> If the initial guess is known to be close to the actual solution is
% best to put that value to about 5. However, if its too far and the
% maximum possible value for the derivative it's too low the convergence is
% slower.
%
% 3 -> Should be an option. The algorithm tries to use low order derivatives
% because they're less expensensive, however, sometimes, it fails to
% converge, so it uses higher order derivatives, that are more precise.