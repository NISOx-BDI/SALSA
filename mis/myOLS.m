function [Bhat,Yhat,res,se] = myOLS(Y,X)

    Y = Y(:); % make it a column 
    T = numel(Y); 

    if size(X,1) ~= T; X = X'; end

    Bhat  = X\Y; % \Beta = X^{+}Y
    Yhat  = X*Bhat; %
    res   = Y-Yhat; % sanitycheck: d.resid(1:5)
    COVX  = (X'*X);
    df    = T-rank(X);
    s2    = (res'*res)./df;
    se    = sqrt(s2./diag(COVX));
end