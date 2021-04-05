function s=newRaps(fun,gradfun,s)
% Finds s when F(s)=0 by Newton Rapshon
% All input are function handles except for the last one. 
% Inputs: fun           =      F(s)
%              gradfun    =     gradient of F(s)
%              s                =      initial gues of F(s)=0
% Outputs: s

% Initial old value
sOld = inf;

% Check if user suplied a function
if isempty(gradfun)
    gradAnalitic=false; 
else
    gradAnalitic=true;
end

% Use finite difference grad
if ~gradAnalitic
    delt=1e-9; % should be an option
    f=fun;
    gradfun=@(s)(f(s-2*delt)-8*f(s-delt)+8*f(s+delt)-f(s+2*delt))/(12*delt);
end

% Main loop
while abs(sOld-s) >= 1e-9 % should be an option
    sOld = s;
    s = nraps(fun,gradfun,s);
end
end

function s=nraps(f,gradfun,s)
    feval=f(s);
    gradf=gradfun(s);
    s=s-feval/gradf;
end