function L = brownian_bridge_bound(epsilon,variance_rate,interval_length,y_difference)
%{ 
tiny convenience function to give a number L such that a Brownian bridge
with variance rate variance_rate starting at (0,0) and ending at
(interval_length,start_end_difference) will exceed L with probability
epsilon. This also works if the bridge starts at (0,start_end_difference)
and ends at (interval_length,0).
%}
a = y_difference;
t = interval_length*variance_rate;
L = (a+sqrt(a^2-2*t*log(epsilon)))/2;
end