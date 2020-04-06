function [cBhat,Yhat,res,stat] = myOLS(Y,X,contrast)

    % scalar: lower-case
    % vector/matrix: capital
    
    Y = Y(:); % make it a column 
    t = numel(Y); 

    if size(X,1) ~= t; X = X'; end

    p           = rank(X);
    stat.Bhat   = X\Y; % \Beta = X^{+}Y
    Yhat        = X*stat.Bhat; %
    res         = Y-Yhat; % sanitycheck: d.resid(1:5)
    COVX        = (X'*X);
    stat.df     = t-p;
    
    
    if ~exist('contrast','var')
        contrast    = zeros(1,p);
        contrast(2) = 1;
    end
        
    stat.mse   = sum(res.^2)/(t-p); 
    stat.se    = sqrt(stat.mse*(contrast*pinv(COVX)*contrast'));
    
    cBhat      = contrast*stat.Bhat;
    stat.tval  = cBhat./stat.se;
    stat.pval  = 2*myOLS_tcdf(-abs(stat.tval),stat.df);
    
    rss = sum(res.^2);
    % these two don't match what Matlab's fitglm.m outputs
    stat.AIC   = 2*(p+1)+t*log(rss/t);
    stat.BIC   = log(t)*(p+1)+t*log(sum(res.^2)/t);
    
end


function p = myOLS_tcdf(x,n);
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

% WOW~
% beautiful for teaching: 
% http://www.karenkopecky.net/Teaching/eco613614/Matlab%20Resources/OLS.pdf

