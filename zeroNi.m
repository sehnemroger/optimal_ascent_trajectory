function x0=zeroNi(f,a,b)
% Find initial estimative of x for f(x)=0 assuming a<x<b 
% and assuming multiplicity m=2.
% This algorithm works with a non interative method. See refs.
% Inputs: f -> Function handle of f
%              a -> Left bound on x
%              b -> Right bound on x
% Output: x0 -> zero aprroximation. 

% Based on "Transformation methods for finding multiple roots of nonlinear
% equations" by Beong In Yun. 
%
% MODIFICATION TO THE METHOD
% fe is now the median derivative far from the zero and is the derivative close to it.
%
% Implemented and modified by Róger M. Sehnem, 2020.
fe=@(x) feps(f,x);
x=-5:.001:5; for i=1:length(x), fep(i)=fe(x(i)); end; plot(x,fep);
Hx=@(x) hyperTangent(f,fe,x);
for i=1:length(x), Hxp(i)=Hx(x(i)); end; plot(x,Hxp);
Ihf=intTrapz(Hx,a,b);
x0=0.5*(a+b)-Ihf/2;
end

function fe=feps(f,x)
    % transformation function that make f have a simple zero
    % and allow local minima.
    % fe tends to the derivative when f(x)->0 and tends to the median
    % gradient with deltap=c when f(x)-> inf or -inf, giving some average
    % derivative.s
    deltax=1e-9; feval=f(x); c=.1; % size of deltap far from the zero, should be an option.
    deltap=deltax+((c-deltax)*feval^2/(1+feval^2));
    fe=(f(x+deltap)-feval)/deltap;
end

function Hx=hyperTangent(f,fe,x)

if f(x) ~= 0
    Hx=tanh(30*fe(x));
else
    Hx=0;
end
end

function int=intTrapz(f,a,b)
% Integrate with trapezoidal rule 

IDv=50; % Interval divisions, Should be an option
x=linspace(a,b,IDv);
h=x(2)-x(1); int=0; feval=f(x(1));
for i=2:length(x)
    fevalm1=feval;
    feval=f(x(i));
    int=int+h*(fevalm1+feval)/2;
end

end