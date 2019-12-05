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
function arrival_list = next_arrival_ihppG_delayV2(params)
% internal parameter for wiener process
epsilon = 1E-6;
% params
c= params.c;
t0 = params.t0;
num_arrivals = params.num_arrivals;
latest_arrival = params.latest_arrival;
d = params.d;
if(isfield(params,'wiener'))
    wiener = params.wiener;% if true, add brownian bridges on to the hash rate function. This also means that the format of lambda is a struct with lambda.mat and lambda.blockrate_f matrix with lambda(1,:) being times and lambda(2,:) being log(hashrate) values, instead of a function
    if(wiener)
        var_rate = params.var_rate;
    end
    ar_piece_x = params.ar_piece_x;
    ar_piece_y = params.ar_piece_y;
else
    wiener = 0;
    subinterval_length = params.subinterval_length;
    lambda = params.lambda;
    bound = params.bound;
end

t1 = t0;
arrival_list = zeros(1,num_arrivals);
num_accepted = 0;
while(num_accepted<num_arrivals)
    if(wiener)% wiener process uses the subintervals defined by ar_piece_x
        if(isempty(ar_piece_x))
            error('ar_piece_x is empty');
        end
        if(t1<ar_piece_x(1)) % if outside the boundary of the domain, extend the hashrate value into the 
            t2 = ar_piece_x(1);
            gB1 = exp(ar_piece_y(1));
            gB2 = exp(ar_piece_y(1));
        elseif(t1>=ar_piece_x(end))
            t2 = t1+ar_piece_x(end)-ar_piece_x(end-1);
            gB1 = exp(ar_piece_y(end));
            gB2 = exp(ar_piece_y(end));
        else
            start_index = find(ar_piece_x==t0);
            if(isempty(start_index))
                error('times not aligned with ar_piece_x');
            end
            t2 = ar_piece_x(start_index+1);
            gB1 = exp(ar_piece_y(start_index));
            gB2 = exp(ar_piece_y(start_index+1));
        end
        subinterval_length = t2-t1;
        current_arrivals = geometric_brownian_poisson(t1,t2,gB1,gB2,var_rate,epsilon);
    else
        current_arrivals = next_arrival_ihppGV(t1,lambda,bound,subinterval_length);% next set of arrivals to test and add
    end
    rejected = zeros(size(current_arrivals));
    for i = 1:numel(current_arrivals)
        rejected(i) = rand>1-exp(-c*(current_arrivals(i)-latest_arrival)-d);
        if(rejected(i)==0)
            latest_arrival = current_arrivals(i);
        end
        %fprintf(1,'rej\n');
    end
    current_arrivals = current_arrivals(rejected==0);
    curr_num_accepted = min(sum(rejected==0),num_arrivals-num_accepted);
    arrival_list((num_accepted+1):(num_accepted+curr_num_accepted)) = current_arrivals(1:curr_num_accepted);
    num_accepted = num_accepted+curr_num_accepted;
    t1 = t1+subinterval_length;
end
end