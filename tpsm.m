function xn=tpsm(func,x)
% Finds s when F(s)=0 based on the Modified Three-Point Secant Method.
%
% Inputs: fun           =      F(s) {function handle}
%              x                =      initial gues of F(s)=0
% Outputs: s


% Method by Ababu T. Tiruneh in "A modified three-point Secant method
%with improved rate and characteristics of convergence".
%
% Implemented by Róger M. Sehnem, 2020.
%% Main loop
counter=0; % initialize the counter
max_error=1e-9; % should be an option
actual_error=inf;
while  actual_error >= max_error
    if counter < 1
        feval=func(x);
        xo1=x-.1;
        xo2=xo1-.1;
        fevalm1=func(xo1);
        fevalm2=func(xo2);
    end
    xn = method(x,feval,fevalm1,fevalm2,xo1,xo2);
    counter=counter+1;
    fevalm2=fevalm1;
    fevalm1=feval;
    xo2=xo1;
    xo1=x;
    x=xn;
    feval=func(x);
    actual_error=abs(feval);
end

counter
end

function xn=method(x,feval,fevalm1,fevalm2,xo1,xo2)

xn=xo2-fevalm2*(feval-fevalm1)/((feval-fevalm2)/(x-xo2)*...
    (feval-fevalm1)-feval*((feval-fevalm2)/(x-xo2)-(fevalm1-fevalm2)/(xo1-xo2)));

end

