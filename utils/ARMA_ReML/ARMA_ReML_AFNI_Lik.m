function LL = ARMA_ReML_AFNI_Lik(a,b,Y,X)
%tic
    m   = size(X,2);
    n   = size(Y,1);
    %pk  = @(a,b,k) (((b+a)*(1+a*b))./(1+2*a*b+b.^2)).*(((b+a)*(1+a*b).*a.^(k-1))./(1+2*a*b+b.^2));
    %k   = 1:n;
    %ppk = pk(a,b,k-1);
    
    %ppk(abs(ppk)<0.01) = 0; 
    %R   = toeplitz([ppk]);
%   a,b
%   R(1:5,1:5)
        
    %pk0 = ((1+ b^2 + 2*a*b)   ./ (1-a^2));
    %pk  = (((b+a) .* (1+a*b)) ./ (1-a^2));
    %pk  = ((a+b)*(1+a*b))./(1+2*a*b+b^2);
    
% Toeplitz
    %A    = [pk pk0 pk];
    %A    = A./pk0;
    %R    = full(spdiags(ones(n,1)*A,(-1:1),n,n));
    %R = ARMACovMat([a,b],n,1,1);
    %R = R./R(1,1); 
    
%      C     = chol(R);
%      Cnegt = inv(C');
%      %D     = Cnegt*X; %size(D)
%      RR = qr(Cnegt*X);
%      RR(1:5,:)
%    
%     %P = Cnegt-Cnegt*X*
%     
%      2*log(det(D))

%     logdetXtRX = 2*sum(log(diag(D)));
%     logdetR    = 2*sum(log(diag(C)));

    %% ---- works fine
    R = ARMACovMat([a,b],n,1,1);
    R(abs(R)<0.01) = 0; 
    %R(1:5,1:5)
    
    invR            = R\eye(n);
    Xt              = X';
    Yt              = Y';
    XtinvR          = Xt*invR; %Xt/R; %Xt*invR
    XtinvRX         = XtinvR*X;
    XinvXtinvRX     = X/XtinvRX; %X*inv(XtinvRX)
    XtinvR          = Xt*invR; %Xt/R; %Xt*invR
    invRXinvXtinvRX = R\XinvXtinvRX; %invR*XinvXtinvRX
    
    P    = invR - invRXinvXtinvRX*XtinvR;
    YtPY = (Yt*P*Y);

    logYtPY    = log(YtPY);
    logdetR    = log(det(R));
    logdetXtRX = log(det(XtinvRX));
    LL   = (n-m).*logYtPY+logdetR+logdetXtRX-log(det(Xt*X));

end