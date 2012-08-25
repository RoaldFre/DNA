% Copyright (C) 1992-1994 Richard Shrager
% Copyright (C) 1992-1994 Arthur Jutan
% Copyright (C) 1992-1994 Ray Muzic
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; If not, see <http://www.gnu.org/licenses/>.

function prt=dfdp(x,f,p,dp,func,bounds)
% function prt = dfdp (x, f, p, dp, func[, bounds])
% numerical partial derivatives (Jacobian) df/dp for use with leasqr
% --------INPUT VARIABLES---------
% x=vec or matrix of indep var(used as arg to func) x=[x0 x1 ....]
% f=func(x,p) vector initialsed by user before each call to dfdp
% p= vec of current parameter values
% dp= fractional increment of p for numerical derivatives
%      dp(j)>0 central differences calculated
%      dp(j)<0 one sided differences calculated
%      dp(j)=0 sets corresponding partials to zero; i.e. holds p(j) fixed
% func=string naming the function (.m) file
%      e.g. to calc Jacobian for function expsum prt=dfdp(x,f,p,dp,'expsum')
% bounds=two-column-matrix of lower and upper bounds for parameters
%      If no 'bounds' options is specified to leasqr, it will call
%      dfdp without the 'bounds' argument.
%----------OUTPUT VARIABLES-------
% prt= Jacobian Matrix prt(i,j)=df(i)/dp(j)
%================================

m=size(x,1); if (m==1), m=size(x,2); end  %# PAK: in case #cols > #rows
n=length(p);      %dimensions
if (nargin < 6)
  bounds = ones (n, 2);
  bounds(:, 1) = -Inf;
  bounds(:, 2) = Inf;
end
prt=zeros(m,n);       % initialise Jacobian to Zero
del = dp .* p; %cal delx=fract(dp)*param value(p)
idx = p == 0;
del(idx) = dp(idx); %if param=0 delx=fraction
idx = dp > 0;
del(idx) = abs (del(idx)); % not for one-sided intervals, changed
				% direction of intervals could change
				% behavior of optimization without bounds
min_del = min (abs (del), bounds(:, 2) - bounds(:, 1));
  for j=1:n
    ps = p;
    if (dp(j)~=0)
      if (dp(j) < 0)
	ps(j) = p(j) + del(j);
	if (ps(j) < bounds(j, 1) || ps(j) > bounds(j, 2))
	  t_del1 = max (bounds(j, 1) - p(j), - abs (del(j))); %
				%non-positive
	  t_del2 = min (bounds(j, 2) - p(j), abs (del(j))); %
				%non-negative
	  if (- t_del1 > t_del2)
	    del(j) = t_del1;
	  else
	    del(j) = t_del2;
	  end
	  ps(j) = p(j) + del(j);
	end
	prt(:, j) = (feval (func, x, ps) - f) / del(j);
      else
	if (p(j) - del(j) < bounds(j, 1))
	  tp = bounds(j, 1);
	  ps(j) = tp + min_del(j);
	elseif (p(j) + del(j) > bounds(j, 2))
	  ps(j) = bounds(j, 2);
	  tp = ps(j) - min_del(j);
	else
	  ps(j) = p(j) + del(j);
	  tp = p(j) - del(j);
	  min_del(j) = 2 * del(j);
	end
	f1 = feval (func, x, ps);
	ps(j) = tp;
        prt(:, j) = (f1 - feval (func, x, ps)) / min_del(j);
      end
    end
  end
