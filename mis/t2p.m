
function [pval,pID,pN] = t2p(tval,df,alev)

if nargin<3
    alev = 0.05; 
end 

ltp   = myOLS_tcdf(tval,df); % tcdf is not available on bmrc Ovtave 
%zval  = sqrt(2)*erfinv(2*ltp-1); % does what norminv do; norminv is not available on bmrc Octave
pval  = 1-ltp;

[pID,pN] = myFDR(pval,alev);

end

function [pID,pN] = myFDR(p,alpha)
% 
% p   - vector of p-values
% alpha   - False Discovery Rate level
%
% pID - p-value threshold based on independence or positive dependence
% pN  - "Nonparametric" (any covariance structure) p-value threshold
%______________________________________________________________________________
% Based on FDR.m     1.4 Tom Nichols 02/07/02

    p = sort(p(:));
    V = length(p);
    I = (1:V)';

    cVID = 1;
    cVN = sum(1./(1:V));

    pID = p(max(find(p<=I/V*alpha/cVID)));
    pN  = p(max(find(p<=I/V*alpha/cVN)));

end


function p = myOLS_tcdf(x,n)
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

end