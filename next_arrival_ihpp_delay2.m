% Calculates the next arrival in an ihpp modified by introducing an 
% exponentially reduced arrival rate for a period after an arrival.
% Uses rejection method.
% c is delay parameter, arrival rate is reduced by a factor of 1-exp(-ct) t
% seconds after an arrival
% lambda = exp(at+b)
% t0 is the time of the most recent arrival
function t1 = next_arrival_ihpp_delay2(params)
c = params.c;
a = params.a;
b = params.b;
t0 = params.t0;
if(isfield(params,'latest_arrival'))
    latest_arrival = params.latest_arrival;
else
    latest_arrival = t0;
end
if(isfield(params,'d'))
    d = params.d;
else
    d = 0;
end
rejected = 1;
t1 = t0;
while(rejected)
    t1 = next_arrival_ihpp(a,b,t1);
    rejected = rand>1-exp(-c*(t1-latest_arrival)-d);
end
end