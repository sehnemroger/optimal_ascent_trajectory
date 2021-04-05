function xn=nsm(func,grad,x)
% Finds s when F(s)=0 based on the Newton-Secant Method. Use analytical
% gradients for improved convergence rate.
%
% Inputs: fun           =      F(s) {function handle}
%              grad          =     dF(s)/ds {function gandle}
%              x                =      initial gues of F(s)=0
% Outputs: s
%
% When using analytical derivatives the algorithm uses a maximum value for
% them, this improve convergence. Source: by experience.
%
% When using the algorithm without the analytical derivatives, there is
% another "option", cMax, that controls how long should the algorithm use
% low order (less expensive) derivatives.
%
% Those options can be easily implemented, but, by now, they're hard
% coded options.

% Method by Gustavo Fernández-Torres in "A Novel Geometric Modification 
% to the Newton-Secant Method to Achieve Convergence of Order 
% 1 + ?2 and Its Dynamics".
%
% Implemented by Róger M. Sehnem, 2020.

%% Decides if it's using analytical derivatives or not
if isempty(grad)
    analytic=false;
else
    analytic=true;
end
%% Main loop
counter=0; % initialize the counter
max_error=1e-9; % should be an option
actual_error=inf;
if analytic
    while  actual_error >= max_error
        if counter < 1
            xo=x;
            fevalm1=func(xo);
            gradf=grad(xo);
            x=xo-fevalm1/gradf;
            feval=func(x);
        end
        if actual_error >= 1e-3 && abs(x-xo) <= 1e-3
            x=xo-rand;
        end            
        xn = nsmAnalytic(x,grad,feval,fevalm1,xo);
        counter=counter+1;
        fevalm1=feval;
        xo=x;
        x=xn;
        feval=func(x);
        actual_error=abs(feval);
%         if counter >= 100
%             break;
%         end
    end
else
    while  actual_error >= max_error
        if counter < 1
            feval=func(x);
            xo=x-1e-3;
            fevalm1=func(xo);
        end
        xn = nsmFiniteDiff(x,feval,fevalm1,xo,func,counter);
        counter=counter+1;
        fevalm1=feval;
        xo=x;
        x=xn;
        feval=func(x);
        actual_error=abs(feval);
    end
end
% counter
end

function xn=nsmAnalytic(x,grad,feval,fevalm1,xo)
%% Calculate function values and derivatives

gradf=grad(x);
% set maximum absolute value for derivatives
% gradmax=50; % -> 2
% gradmin=-gradmax;
% gradf=min(max(gradf,gradmin),gradmax);
%

xn=x-feval*(feval-fevalm1)*(x-xo)/(feval*(feval-fevalm1)-fevalm1*gradf*(x-xo));
end

function xn=nsmFiniteDiff(x,feval,fevalm1,xo,f,counter)
% Function variables are have not  the same name because i'm too lazy
delt=.5e-8;
cMax=500; % -> 3
%
if counter <= cMax
    % Low order estimate
    fma1=f(x+1*delt); %fma2=f(x+2*delt);
    gradf=(fma1-feval)/(delt);
%     gradf2=(feval-2*fma1+fma2)/(delt^2);
else
    % High order estimate
    fme2=f(x-2*delt); fme1=f(x-1*delt); fma1=f(x+1*delt); fma2=f(x+2*delt);
    gradf=(fme2-8*fme1+8*fma1-fma2)/(12*delt);
end
%
xn=x-feval*(feval-fevalm1)*(x-xo)/(feval*(feval-fevalm1)-fevalm1*gradf*(x-xo));
end

