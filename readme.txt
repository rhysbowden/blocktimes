A project for simulating block mining times in the Bitcoin blockchain based on the paper "Modelling and analysis of block arrival times in the Bitcoin blockchain" by Bowden et al.

The main file is blocktime_sim_delay.m. Here is an explanation of the input parameters.

[block_times,difficulty_history,varargout] = blocktime_sim_delay(num_blocks,hashrate_params,starting_difficulty,starting_time,starting_segment_position,latest_epoch_time,latest_arrival)
Simulate from a model for the blockchain

num_blocks is the number of blocks to calculate arrival times for (not
  including the zeroth block, at time 0)
hashrate_params are specific to the type of model, see below.
  hashrate_params.c and hashrate_params.d are the parameters for delay,
  see the paper Bowden et al., "Modelling block arrivals in the Bitcoin
  blockchain". It is not recommended using 
starting_difficulty is the average number of hashes that must be checked
  to mine a block successfully (2^256/target), rather than the common 
  definition of difficulty which is approximately 2^224/target. The 
  higher the difficulty, the harder it is to mine blocks.
starting time is the start time of the simulation.
starting_segment_position is how many blocks there have been since the
  latest difficulty change.
latest_epoch_time is the most recent time at which a difficulty change
  occurred.
latest_arrival is the most recent time at which a block arrived.

hashrate_params.hashrate_type determines which type of model is 
  simulated, the default is 'exp'
If hashrate_params.hashrate_type='exp' then the hashrate H(t) = exp(at+b)
  where a is given by hashrate_params.a and b is given by
  hashrate_params.b.
If hashrate_params.hashrate_type='deterministic exp' the model is the 
  same as 'exp' but difficulty adjustments are deterministically 
  calculated ahead of time, and the whole process is a non-homogeneous 
  Poisson process.
If hashrate_params.hashrate_type='general' then the hashrate is given by
  a function handle hashrate_params.hashrate_f(t). bound(t1,t2) is the 
  maximum value of the hashrate function in the interval [t1,t2]. If the 
  hashrate function is increasing then it will just be hashrate_f(t2).
  This is required to use the rejection method to sample block times.
  It is recommended that you used 'generalV' instead of 'general' since
  it is typically faster.
If hashrate_params.hashrate_type='generalV', then the model is the same
  as in 'general'. hashrate_params.subinterval_length is the length of
  each period of sampling times. The default is 1 day. 
  However, if also hashrate_params.wiener==1, then the
  hashrate will be specified by hashrate_params.hr_piece_x and 
  hashrate_params.hr_piece_y which give the hashrate exp(hr_piece_y) at 
  the times given in hr_piece_x, and the logarithm of the hashrate 
  function is interpolated between these times with a Brownian bridge. 
  The variance rate of these Brownian bridges is given by 
  hashrate_params.var_rate. hashrate_params.gB_curr is the hashrate at
  the starting time of the simulation.
