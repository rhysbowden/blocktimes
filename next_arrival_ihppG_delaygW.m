% Calculates the next num_arrivals arrivals in an ihpp modified by introducing an
% exponentially reduced arrival rate for a period after an arrival.
% Uses rejection method.
% c is delay parameter, arrival rate is reduced by a factor of 1-exp(-ct) t
% seconds after an arrival
% lambda() is a nonparametric function bounded above on (x,y) by bound(x,y)
% t0 is the time of the most recent arrival
% based on next_arrival_ihppG_delay but modified to do multiple arrivals at
% once
% version 2 uses struct as input, and can take a value of d, so that
% arrival rate is instead reduced by a factor of (1-exp(-ct-d)) t seconds
% after arrival.
% 8 Apr 2019: added wiener process version
% version gW: geometric brownian motion, needs another parameter, the
% starting hash rate.
function [arrival_list,t2,gB2] = next_arrival_ihppG_delaygW(params,gB0)
% internal parameter for wiener process
epsilon = 1E-6;
% params
c= params.c;
t0 = params.t0;
num_arrivals = params.num_arrivals;
latest_arrival = params.latest_arrival;
d = params.d;
var_rate = params.var_rate;
ar_piece_x = params.ar_piece_x;
ar_piece_y = params.ar_piece_y;% log of arrival rate
if(isempty(ar_piece_x))
    error('ar_piece_x is empty');
end

t1 = t0;% running start time
gB1 = gB0;
arrival_list = zeros(1,num_arrivals);
num_accepted = 0;
while(num_accepted<num_arrivals)
    % -------- Poisson arrivals part -------------------------------
    if(t1>=ar_piece_x(end))
        t2 = t1+ar_piece_x(end)-ar_piece_x(end-1);
        gB2 = exp(ar_piece_y(end));
    else
        t2_ind = find(ar_piece_x>t1,1,'first');
        t2 = ar_piece_x(t2_ind);
        gB2 = exp(ar_piece_y(t2_ind));
    end
    [current_arrivals,arrival_rates] = geometric_brownian_poisson(t1,t2,gB1,gB2,var_rate,epsilon);% arrival rates is the vector of arrival rates, one for each point in current_arrivals.
    
    % -------- Rejection due to propagation delay part.--------------
    rejected = zeros(size(current_arrivals));
    for i = 1:numel(current_arrivals)
        rejected(i) = rand>1-exp(-c*(current_arrivals(i)-latest_arrival)-d);
        if(rejected(i)==0)
            latest_arrival = current_arrivals(i);
        end
        %fprintf(1,'rej\n');
    end
    current_arrivals = current_arrivals(rejected==0);
    arrival_rates = arrival_rates(rejected==0);
    curr_num_accepted = min(sum(rejected==0),num_arrivals-num_accepted);
    if(sum(rejected==0)>(num_arrivals-num_accepted))% if we are cutting off the end of (t1,t2) due to hitting the max number of arrivals
       t2 = current_arrivals(curr_num_accepted);
       gB2 = arrival_rates(curr_num_accepted);
    end
    arrival_list((num_accepted+1):(num_accepted+curr_num_accepted)) = current_arrivals(1:curr_num_accepted);
    num_accepted = num_accepted+curr_num_accepted;
    t1 = t2;
    gB1 = gB2;
end
end