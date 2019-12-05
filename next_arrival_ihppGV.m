% next_arrival_ihppG
% rejection method for calculating arrivals in next period for an IHPP with
% intensity function lambda(t)
% runs from t0 to t0+subinterval_length
% bound is a function that takes the ends of an interval [x,y] and is
%   at least lambda(t) for all t in [x,y] 
function arrivals = next_arrival_ihppGV(t0,lambda,bound,subinterval_length)
vectorised = 1;
tend = t0+subinterval_length;% end of interval
L = bound(t0,tend);%upper limit of lambda over the interval
num_arrivals = random('Poisson',bound(t0,tend)*(tend-t0));
arrivals = random('Uniform',t0,tend,num_arrivals,1);
arrivals = sort(arrivals,'ascend');
if(vectorised)% if lambda accepts vectors as input
    arrivals = arrivals(rand(num_arrivals,1)<(lambda(arrivals)/L));
    if(any(lambda(arrivals)/L>1))
        error('bound(x,y) is not always greater than lambda');
    end
else
    accepted = zeros(numel(arrivals),1);
    for i =1:numel(arrivals) % loop needed because lambda isn't vectorised
        accepted(i) = rand<(lambda(arrivals(i))/L);
    end
    arrivals=arrivals(accepted>0);
end
end