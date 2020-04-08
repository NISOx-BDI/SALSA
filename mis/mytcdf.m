function p = mytcdf(x,n);
% TCDF returns student cumulative distribtion function
%
% cdf = tcdf(x,DF);
%

% check size of arguments
if all(size(x)==1)
	x = repmat(x,size(n));
elseif all(size(n)==1)
	n = repmat(n,size(x));
elseif all(size(x)==size(n))
	;	%% OK, do nothing
else
	error('size of input arguments must be equal or scalar')
end; 	
		
% allocate memory
p = zeros(size(x));
p((x==Inf) & (n>0)) = 1;

% workaround for invalid arguments in BETAINC
ix   = isnan(x) | ~(n>0);
p(ix)= NaN;

ix    = (x > -Inf) & (x < Inf) & (n > 0);
p(ix) = betainc (n(ix) ./ (n(ix) + x(ix).^2), n(ix)/2, 1/2) / 2;

ix    = find(x>0);
p(ix) = 1 - p(ix);

% shape output
p = reshape(p,size(x));
