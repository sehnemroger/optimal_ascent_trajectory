function x0=iAproxm2(f,a,b)
% from "Transformation methods for finding multiple roots of nonlinear
% equations" by Beong In Yun.
%
% Implemented by Róger M. Sehnem, 2020.

epsl=1e-9;
fe=@(x) f(x+epsl*f(x))-f(x);
Hx=@(x) hyperTangent(f,fe,x);
Ihf=intQuad(Hx,a,b);
x0=0.5*(a+b)-Ihf/2;
end

function Hx=hyperTangent(f,fe,x)

if f(x) ~= 0
    Hx=tanh(1./fe(x));
else
    Hx=0;
end
end

function int=intQuad(f,a,b)
% Integrate with trapezoidal rule 

IDv=50; % Interval divisions
x=linspace(a,b,IDv);
h=x(2)-x(1); int=0; feval=f(x(1));
for i=2:length(x)
    fevalm1=feval;
    feval=f(x(i));
    int=int+h*(fevalm1+feval)/2;
end

end