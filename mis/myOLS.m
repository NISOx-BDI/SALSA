function [cbhat,Yhat,RES,stat] = myOLS(Y,X,contrast)
    % intercept should always be the first column
    % scalar: lower-case
    % vector/matrix: capital
    % 
    % SA, Ox, 2020
    %
    verbose = 0;
    
    if size(Y,1)~=size(X,1); error('myOLS:: check the inputs!'); end; 
    if size(X,1)<size(X,2); error('myOLS:: The model is not stable.'); end; 
    
    
    t           = size(X,1);
    p           = rank(X);
    
    if verbose 
        disp(['voxels: ' num2str(size(Y,2)) ', timepoints: ' num2str(t) ', rankX: ' num2str(p) ])
    end
    
    stat.Bhat   = X\Y; % \Beta = X^{+}Y
    Yhat        = X*stat.Bhat; %
    RES         = Y-Yhat; % sanitycheck: d.resid(1:5)
    COVX        = (X'*X);
    stat.df     = t-p;
    
    if ~exist('contrast','var') % assuming the first column is the intercept
        contrast    = zeros(1,p);
        contrast(2) = 1;
    end
    
    VarYhat    = sum(Yhat.^2);
    VarY       = sum(Y.^2);
    R2         = VarYhat./VarY; % R-squared
    stat.f     = (R2*(t-p))./((1-R2).*(p-1)); %F statsitics
    stat.fp    = 1-fcdf(stat.f,(p-1),(t-p));
    
    
    stat.mse   = sum(RES.^2)./(t-p);
    stat.se    = sqrt(stat.mse*(contrast*pinv(COVX)*contrast'));
    
    cbhat      = contrast*stat.Bhat;
    stat.tval  = cbhat./stat.se;
    stat.ltp   = myOLS_tcdf(stat.tval,stat.df); % tcdf is not available on bmrc Ovtave 
    stat.zval  = sqrt(2)*erfinv(2*stat.ltp-1); % does what norminv do; norminv is not available on bmrc Octave
    stat.pval  = 1-stat.ltp;
        
    
    % these two don't match what Matlab's fitglm.m outputs
%     rss        = sum(RES.^2);
%     stat.AIC   = 2*(p+1)+t*log(rss/t);
%     stat.BIC   = log(t)*(p+1)+t*log(sum(RES.^2)/t);
    
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

% ======= SPM ==================
% Contrast
% C = [1 0];
% % t statistic and significance test
% RSS   = sum((Y - X*B).^2);
% MRSS  = RSS / df;
% SE    = sqrt(MRSS*(C*pinv(X'*X)*C'));
% t     = C*B./SE;
% ltp   = spm_Tcdf(t, df); % lower tail p
% Z     = spm_invNcdf(ltp);
% p     = 1-ltp;           % upper tail p 
% ==============================

% WOW~
% beautiful for teaching: 
% http://www.karenkopecky.net/Teaching/eco613614/Matlab%20Resources/OLS.pdf

