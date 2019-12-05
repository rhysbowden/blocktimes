% Calculates the next arrival in an ihpp modified by introducing an 
% exponentially reduced arrival rate for a period after an arrival.
% Uses rejection method.
% c is delay parameter, arrival rate is reduced by a factor of 1-exp(-ct) t
% seconds after an arrival.
% lambda = exp(at+b)
% t0 is the time of the most recent arrival
function t1 = next_arrival_ihpp_delay(c,a,b,t0,varargin)
if(nargin==4)
    latest_arrival = t0;
elseif(nargin==5)
    latest_arrival = varargin{1};
else
    error('Wrong number of arguments for next_arrival_ihpp_delay');
end
rejected = 1;
t1 = t0;
while(rejected)
t1 = next_arrival_ihpp(a,b,t1);
rejected = rand>1-exp(-c*(t1-latest_arrival));
%fprintf(1,'rej\n');
end
end