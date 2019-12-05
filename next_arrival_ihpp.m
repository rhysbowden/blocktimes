% inverse transform method for calculating next arrival in an IHPP
% with lambda(t) = exp(at+b)
function t1 = next_arrival_ihpp(a,b,varargin)
if(nargin>2)
    t0 = varargin{1};
else
    t0 = 0;
end
if(a~=0)
    t1 = log(exp(a*t0)-a*log(rand)/exp(b))/a;
else
    t1 = random('exp',1/exp(b))+t0;
end
end