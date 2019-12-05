% [block_times,difficulty_history,varargout] = blocktime_sim_delay(num_blocks,hashrate_params,starting_difficulty,starting_time,starting_segment_position,latest_epoch_time,latest_arrival)
% Simulate from a model for the blockchain
%
% num_blocks is the number of blocks to calculate arrival times for (not
%   including the zeroth block, at time 0)
% hashrate_params are specific to the type of model, see below.
%   hashrate_params.c and hashrate_params.d are the parameters for delay,
%   see the paper Bowden et al., "Modelling block arrivals in the Bitcoin
%   blockchain". It is not recommended using 
% starting_difficulty is the average number of hashes that must be checked
%   to mine a block successfully (2^256/target), rather than the common 
%   definition of difficulty which is approximately 2^224/target. The 
%   higher the difficulty, the harder it is to mine blocks.
% starting time is the start time of the simulation.
% starting_segment_position is how many blocks there have been since the
%   latest difficulty change.
% latest_epoch_time is the most recent time at which a difficulty change
%   occurred.
% latest_arrival is the most recent time at which a block arrived.
%
% hashrate_params.hashrate_type determines which type of model is 
%   simulated, the default is 'exp'
% If hashrate_params.hashrate_type='exp' then the hashrate H(t) = exp(at+b)
%   where a is given by hashrate_params.a and b is given by
%   hashrate_params.b.
% If hashrate_params.hashrate_type='deterministic exp' the model is the 
%   same as 'exp' but difficulty adjustments are deterministically 
%   calculated ahead of time, and the whole process is a non-homogeneous 
%   Poisson process.
% If hashrate_params.hashrate_type='general' then the hashrate is given by
%   a function handle hashrate_params.hashrate_f(t). bound(t1,t2) is the 
%   maximum value of the hashrate function in the interval [t1,t2]. If the 
%   hashrate function is increasing then it will just be hashrate_f(t2).
%   This is required to use the rejection method to sample block times.
%   It is recommended that you used 'generalV' instead of 'general' since
%   it is typically faster.
% If hashrate_params.hashrate_type='generalV', then the model is the same
%   as in 'general'. hashrate_params.subinterval_length is the length of
%   each period of sampling times. The default is 1 day. 
%   However, if also hashrate_params.wiener==1, then the
%   hashrate will be specified by hashrate_params.hr_piece_x and 
%   hashrate_params.hr_piece_y which give the hashrate exp(hr_piece_y) at 
%   the times given in hr_piece_x, and the logarithm of the hashrate 
%   function is interpolated between these times with a Brownian bridge. 
%   The variance rate of these Brownian bridges is given by 
%   hashrate_params.var_rate. hashrate_params.gB_curr is the hashrate at
%   the starting time of the simulation.

function [block_times,difficulty_history,varargout] = blocktime_sim_delay(num_blocks,hashrate_params,starting_difficulty,starting_time,starting_segment_position,latest_epoch_time,latest_arrival)
real_difficulty_adjustment = hashrate_params.real_difficulty_adjustment;% 1 if use only last 2015 blocks for difficulty adjustment as per the bug in the real blockchain
% parameters
hashrate_type = hashrate_params.hashrate_type;
dd = 0;
if(isfield(hashrate_params,'d'))
    dd = hashrate_params.d;
end
if(strcmpi(hashrate_type,'exp')) % model 1, exponential random
    hrt = 1;% hashrate type number
    aa = hashrate_params.a;
    bb = hashrate_params.b; % hashrate is exp(at+b)
    cc = hashrate_params.c;
elseif(strcmpi(hashrate_type,'deterministic exp')) % model 2, exponential deterministic version
    hrt=2;
    aa = hashrate_params.a;
    bb = hashrate_params.b; % hashrate is exp(at+b)
    cc = Inf;%hashrate_params.c; % not recommended using c with this model because the difficulty calculation formula is untested.
    duration = hashrate_params.duration; % for deterministic, ignore num_arrivals
elseif(strcmpi(hashrate_type,'general'))% model 3, general random hash rate, individual arrivals (slow)
    hrt = 3;% hashrate type number
    hashrate_f = hashrate_params.hashrate_f;
    bound_f = hashrate_params.bound_f;
    subinterval_length = hashrate_params.subinterval_length;
    cc = hashrate_params.c;
elseif(strcmpi(hashrate_type,'generalV')) % model 3, vectorised version of general hash rate
    hrt = 3.5;% hashrate type number
    hashrate_f = hashrate_params.hashrate_f;
    bound_f = hashrate_params.bound_f;
    cc = hashrate_params.c;
    if(isfield(hashrate_params,'wiener'))
        wiener = hashrate_params.wiener;
    else
        wiener = 0;
    end
    if(wiener)
        hr_piece_x = hashrate_params.hr_piece_x;% grid of times in seconds
        hr_piece_y = hashrate_params.hr_piece_y;% log of hashrate
        var_rate = hashrate_params.var_rate;
        t_curr = starting_time;
        gB_curr = hashrate_params.gB_curr;
        num_arrivals = hashrate_params.num_arrivals;
    else
        if(isfield(hashrate_params,'subinterval_length'))
            subinterval_length = hashrate_params.subinterval_length;
        else
            subinterval_length = 86400; % the default: 1 day in seconds
        end
        num_arrivals = hashrate_params.num_arrivals;
    end
else % assume exp
    hrt = 1;
    aa = hashrate_params.a;
    bb = hashrate_params.b; % hashrate is exp(at+b)
    cc = hashrate_params.c;
end

intended_block_interval = 600;% 10 minutes (in seconds)
intended_segment_time = intended_block_interval*2016;% 2 weeks by default
difficulty_history = zeros(1,num_blocks);% also includes the new difficulty if the simulation ends on a difficulty boundary
block_times = zeros(num_blocks,1);
% use inverse transform method to get inhomogeneous poisson process
% combined with rejection to introduce delay
current_difficulty = starting_difficulty;
current_time = starting_time;
if(hrt==1||hrt==3)% hash rate exponential or general (nonparametric) $ model 1 or 3 (3 is slow this way)
    for i = 1:num_blocks % i is index of next block to arrive
        if(hrt==1)
            if(dd==0) % for backwards compatibility
                block_times(i) = next_arrival_ihpp_delay(cc,aa,bb-log(current_difficulty),current_time,latest_arrival);
            else
                next_arrival_params = struct();
                next_arrival_params.c = cc;
                next_arrival_params.a = aa;
                next_arrival_params.b = bb-log(current_difficulty);
                next_arrival_params.t0 = current_time;
                next_arrival_params.latest_arrival = latest_arrival;
                next_arrival_params.d = dd;
                block_times(i) = next_arrival_ihpp_delay2(next_arrival_params);
            end
        elseif(hrt==3)
            if(dd~=0)
                error('Use vectorised version when d>0');
            end
            block_times(i) = next_arrival_ihppG_delay(cc,current_time,@(x)(hashrate_f(x)/current_difficulty),@(x,y)(bound_f(x,y)/current_difficulty),subinterval_length,latest_arrival);
        end
        latest_arrival = block_times(i);
        current_time = block_times(i);
        difficulty_history(i) = current_difficulty;
        if(mod(i+starting_segment_position,2016)==0)
            if(i<=2016&&~real_difficulty_adjustment)
                % use initial values
                actual_segment_time = block_times(i)-latest_epoch_time;
            else
                % difficulty adjustment
                actual_segment_time = block_times(i)-block_times(i-2016+real_difficulty_adjustment); % 2015 block-intervals version of the real chain
            end
            current_difficulty = current_difficulty*min(max(intended_segment_time/actual_segment_time,0.25),4);
        end
    end
elseif(hrt==2) % deterministic exponential % model 2
    % generate segment boundaries (and difficulties)
    estimated_max_number_of_segments = 1+ceil(duration/(86400*7)); % one more than number of weeks in duration, rounded up
    segment_times = zeros(1,estimated_max_number_of_segments); % do this so we don't have to resize inside loop
    segment_difficulties = zeros(1,estimated_max_number_of_segments);
    i=1;
    current_segment_start = starting_time;
    current_difficulty = starting_difficulty;
    while(current_segment_start<starting_time+duration)
        %next_segment_start = (log((exp(aa*current_segment_start+bb)+cc)*exp(2016*current_difficulty*aa/cc)-cc)-bb)/aa;
        % calculate the end of the current segment/start of next segment by
        %   inverting an integral
        % The indefinite integral is log(cc+exp(aa*tt/bb)/difficulty) and
        % it is evaluated from current_segment_start to next_segment_start.
        % We then rearrange to make next_segment_start the subject of the
        % equation without rounding errors.
        if(cc == Inf)
            %next_segment_start = log(exp(aa*current_segment_start)+2016*current_difficulty*aa/exp(bb))/aa;
            next_segment_start = current_segment_start+log(2016*current_difficulty*aa/(exp(aa*current_segment_start+bb))+1)/aa;
        else
            next_segment_start = (log((exp(2016*aa/cc)*(cc+(1/current_difficulty)*exp(aa*current_segment_start+bb))-cc)*current_difficulty)-bb)/aa;% this might be wrong! RB
        end
        segment_times(i) = next_segment_start;
        segment_difficulties(i) = current_difficulty;
        % update difficulty
        actual_segment_time = next_segment_start - current_segment_start;
        if(real_difficulty_adjustment)
            actual_segment_time = actual_segment_time*2015/2016;
        end
        current_difficulty = current_difficulty*min(max(intended_segment_time/actual_segment_time,0.25),4);
        % update references for next iteration of this loop
        current_segment_start = next_segment_start;
        i=i+1;
        % resize segment_times if necessary
        if(i>estimated_max_number_of_segments)
            segment_times = [segment_times zeros(1,numel(segment_times))]; %#ok<AGROW>
            estimated_max_number_of_segments = estimated_max_number_of_segments*2;
        end
    end
    segment_times = segment_times(1:i-1); % cut off unused space
    segment_difficulties = segment_difficulties(1:i-1);
    num_segments = i-1; % includes partial segment at end (possibly length 0)
    
    % simulate block arrivals
    estimated_max_number_of_blocks = ceil(2016*num_segments*(1+0.1/sqrt(num_segments))); % just a guess for preallocation purposes
    block_times = zeros(1,estimated_max_number_of_blocks);
    j=1; % index of next block
    latest_arrival = starting_time;% the most recent arrival
    t0 = starting_time;% the time up to which we have checked for arrivals
    for i=1:num_segments
        while(1)
            next_arrival = next_arrival_ihpp_delay(cc,aa,bb-log(segment_difficulties(i)),t0,latest_arrival);
            if(next_arrival>segment_times(i)) % if we go over a segment boundary.
                t0 = segment_times(i);
                break;
            else
                if(j>estimated_max_number_of_blocks) % if we go past the end of the block_times array then resize it
                    array_size_increment = ceil(0.1*estimated_max_number_of_blocks);
                    block_times = [block_times zeros(1,array_size_increment)]; %#ok<AGROW>
                    estimated_max_number_of_blocks = estimated_max_number_of_blocks+array_size_increment;
                end
                % record new block time and difficulty
                block_times(j) = next_arrival;
                latest_arrival = next_arrival;
                t0 = next_arrival;
                j=j+1;
                difficulty_history(j) = segment_difficulties(i);
            end
        end
    end
    block_times = block_times(1:(j-1));% check this
elseif(hrt==3.5)% hash rate nonparametric, vectorised arrivals (faster), random difficulty changes % model 3
    i=1;% index of next block to arrive
    while i<=num_blocks
        next_arrival_params = struct;
        next_arrival_params.c = cc;
        next_arrival_params.latest_arrival = latest_arrival;
        next_arrival_params.d = dd;
        arrival_cap = 2016-mod(i-1+starting_segment_position,2016); % number of arrivals still to happen before difficulty change
        if(wiener)
            next_arrival_params.ar_piece_x = hr_piece_x; % t values for the arrival rate function.
            next_arrival_params.ar_piece_y = hr_piece_y - log(current_difficulty);% log of arrival rate
            next_arrival_params.var_rate = var_rate;
            next_arrival_params.num_arrivals = min([num_arrivals,num_blocks-i+1,arrival_cap]);
            next_arrival_params.t0 = t_curr;
            [arrival_list,t_curr,gB_curr] = next_arrival_ihppG_delaygW(next_arrival_params,gB_curr);
        else
            next_arrival_params.t0 = latest_arrival;
            next_arrival_params.num_arrivals = num_arrivals;
            next_arrival_params.subinterval_length = subinterval_length;
            next_arrival_params.lambda = @(x)(hashrate_f(x)/current_difficulty);
            next_arrival_params.bound = @(x,y)(bound_f(x,y)/current_difficulty);
            arrival_list = next_arrival_ihppG_delayV2(next_arrival_params);
        end
        curr_num_arrivals = numel(arrival_list);
        if(curr_num_arrivals>arrival_cap) % if we pass a difficulty change with this set of arrivals
            arrival_list = arrival_list(1:arrival_cap); % then remove all arrivals after the difficulty change
            curr_num_arrivals = arrival_cap;
        end
        block_times(i:i+curr_num_arrivals-1) = arrival_list';
        difficulty_history(i:i+curr_num_arrivals-1) = current_difficulty;
        latest_arrival = arrival_list(end);
        i=i+curr_num_arrivals;
        % difficulty change:
        if(mod(i-1+starting_segment_position,2016)==0)
            if(i<=2017&&~real_difficulty_adjustment)
                % use initial values
                actual_segment_time = block_times(i-1)-latest_epoch_time;
            else
                % difficulty adjustment
                actual_segment_time = block_times(i-1)-block_times(i-2017+real_difficulty_adjustment); % 2015 block-intervals version of the real chain
                % not starting on a block boundary and having
                % real_difficulty_adjustment ==1 will cause an error
            end
            difficulty_change_factor = min(max(intended_segment_time/actual_segment_time,0.25),4);
            current_difficulty = current_difficulty*difficulty_change_factor;
            if(wiener)
               gB_curr = gB_curr*difficulty_change_factor; 
            end
        end
    end
end
if(nargout>2)
    varargout{1} = segment_times;
end
end