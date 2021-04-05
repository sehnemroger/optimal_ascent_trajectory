function xn=secT(func,grad,grad2,x)
% Finds s when F(s)=0 based on the New Secant Type Method. Using analytical
% gradients for improved convergence rate.
%
% Inputs: fun           =      F(s) {function handle}
%              grad          =     dF(s)/ds {function gandle}
%              x                =      initial gues of F(s)=0
% Outputs: s
%

% Method by R. Thukral in "A New Secant-type Method for Solving Nonlinear
% Equations".
%
% Implemented by Róger M. Sehnem, 2020.

%% Main loop
travou=false;
counter=0; % initialize the counter
max_error=1e-9; % should be an option
actual_error=inf;
    while  actual_error >= max_error
        if counter < 1
            xo=x;
            fevalm1=func(xo);
            gradfm1=grad(xo);
            x=xo-fevalm1/gradfm1;
            feval=func(x);
            gradf=grad(x);
        end
%         if actual_error >= 1e-3 && abs(x-xo) <= 1e-5
%             x=x*rand;
%         end
        if counter >= 800
            travou=true;
%             abs(x-xo)
%             x_plot=-15:.01:15; for i=1:length(x_plot), y_plot(i)=func(x_plot(i)); end
%             plot(x_plot,y_plot);
%             hold on
%             plot(x,func(x),'*r');
            break
        end
        xn = secTAnalytic(x,xo,feval,gradf,gradfm1);
        counter=counter+1;
        gradfm1=gradf;
        xo=x;
        x=xn;
        feval=func(x);
        gradf=grad(x);
        actual_error=abs(feval);
    end
    if travou
        xn=step2b(func,grad,grad2,x);
    end
counter
end

function xnew=secTAnalytic(x,xo,feval,gradf,gradfm1)
%% Calculate function values and derivatives
% set maximum absolute value for derivatives
% gradmax=1e4; % -> 2
% gradmin=-gradmax;
% gradf=min(max(gradf,gradmin),gradmax);

xnew=x-2*feval*gradf*(x-xo)/(2*gradf^2*(x-xo)-feval*(gradf-gradfm1));
end