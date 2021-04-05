function s=twoStep(fun,s)
% Finds s when F(s)=0 by the Two Step Iterative Method for Finding Root of a
% Nonlinear Equation.
% Rajesh C. Shah, Rajiv B. Shah
%
% Inputs: fun           =      F(s) {function handle}
%              s                =      initial gues of F(s)=0
% Outputs: s

% Main loop
counter=0; % initialize the counter
while abs(fun(s)) >= 1e-9 % should be an option
    s = twoStepfunc(fun,s,counter);
    counter=counter+1;
end
end

function phi=twoStepfunc(f,phi,counter)
% Two Step Iterative Method for Finding Root of a
% Nonlinear Equation. Rajesh C. Shah, Rajiv B. Shah
delt=.5e-8;
feval=f(phi);
cMax=50;
%
if counter <= cMax
    % Estimativa de baixa ordem
    fma1=f(phi+1*delt); fma2=f(phi+2*delt);
    gradf=(fma1-feval)/(delt); %grad numérico segunda ordem
    gradf2=(feval-2*fma1+fma2)/(delt^2);%grad2 numérico terceira ordem
else
    % Estimativa de derivadas de mais alta ordem
    fme2=f(phi-2*delt); fme1=f(phi-1*delt); fma1=f(phi+1*delt); fma2=f(phi+2*delt);
    gradf=(fme2-8*fme1+8*fma1-fma2)/(12*delt); %grad numérico quinta ordem
    gradf2=-(fme2-16*fme1+30*feval-16*fma1+fma2)/(12*delt^2);%grad2 numérico quinta ordem
end
%
%
m=2; % parametro de escolha
%
% parteRaiz=real(sqrt(1/m^2*feval^2*gradf^2-(1-1/m^2)*feval^3*gradf2));
parteRaiz=sqrt(1/m^2*feval^2*gradf^2-(1-1/m^2)*feval^3*gradf2);
if ~isreal(parteRaiz), parteRaiz=0; end
%                          ou +
yn=phi-(feval*gradf-parteRaiz)/(feval*gradf2+gradf^2);
fevalyn=f(yn);
%
if counter <= cMax
    % Estimativa de baixa ordem
    fma1yn=f(yn+1*delt);
    gradfyn=(fma1yn-fevalyn)/(delt); %grad numérico segunda ordem
    %
else
    % Estimativa de derivadas de mais alta ordem
    fme2yn=f(yn-2*delt); fme1yn=f(yn-1*delt); fma1yn=f(yn+1*delt); fma2yn=f(yn+2*delt);
    gradfyn=(fme2yn-8*fme1yn+8*fma1yn-fma2yn)/(12*delt); %grad numérico quinta ordem
end

phi=yn-fevalyn/gradfyn;
end