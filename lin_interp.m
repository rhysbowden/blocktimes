% linear interpolation with constant values at the ends
function y = lin_interp(x,x_vals,y_vals)
% if(sum(x>x_vals)==0)
%     y = y_vals(1);
% elseif(sum(x<x_vals)==0)
%     y = y_vals(end);
% else
%     ind = find(x_vals<x,1,'last');
%     y = y_vals(ind)+(y_vals(ind+1)-y_vals(ind))*(x-x_vals(ind))/(x_vals(ind+1)-x_vals(ind));
% end
y = interp1(x_vals,y_vals,x);
y(x<x_vals(1)) = y_vals(1);
y(x>x_vals(end)) = y_vals(end);
end