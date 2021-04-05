function [fa2,gfa2,g2fa2]=fa2(x)
% Quadratic penalization for x>0

if x <= 0
    fa2=0;
elseif x>0
    fa2=x^2;
end

if nargout > 1
    if x <= 0
        gfa2=0;
    elseif x > 0
        gfa2=2*x;
    end
    if nargout>2
        if x <= 0
            g2fa2=0;
        elseif x > 0
            g2fa2=2;
        end
    end
end
end