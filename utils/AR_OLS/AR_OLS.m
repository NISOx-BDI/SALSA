function rho = AR_OLS(Y,p)
% Do not turn this into a matrix form. This can be easily translated to
% Python as it is. 
% SA, Ox, 2020

% T=1000;
% model = arima('AR',{0.1},'Constant',0,'Variance',1);
% Y = simulate(model,T);

    if ~exist('p','var') || isempty(p); p=1; end

    T    = numel(Y);
    Ymat = zeros(T-p,p);
    Yb   = zeros(T-p-1,1);
    for c=p:-1:1
        Ymat(:, 1+(p-c)) = Y(1+(p-c):end-c);
    end
        
    Yb  = Y((p+1):end);    
    rho = flip(pinv(Ymat)*Yb);

end