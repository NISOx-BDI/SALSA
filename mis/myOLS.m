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
    stat.pval  = 2*tcdf(-abs(stat.tval),stat.df);
    
    rss = sum(res.^2);
    % these two don't match what Matlab's fitglm.m outputs
    stat.AIC   = 2*(p+1)+t*log(rss/t);
    stat.BIC   = log(t)*(p+1)+t*log(sum(res.^2)/t);
    
end




% WOW~
% beautiful for teaching: 
% http://www.karenkopecky.net/Teaching/eco613614/Matlab%20Resources/OLS.pdf

