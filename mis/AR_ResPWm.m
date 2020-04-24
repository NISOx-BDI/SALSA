function [sqrtmVhalf,posdef] = AR_ResPWm(autocov,AROrd,BiasAdj)
%
% dRESacov : autocovariance 
% Mord : Order required 
% BiasAdj: bias adjustment matrix 
%
% SA, Ox, 2020
%
    if nargin < 3 || isempty(BiasAdj); BiasAdj = eye(AROrd+1); end;

    autocov                 =  autocov(:);
    
    T                       = numel(autocov); % full ACF has one element less than the time series
    dRESac_adj              = (BiasAdj*autocov(1:AROrd+1));
    dRESac_adj              = dRESac_adj./dRESac_adj(1); % make a auto*correlation*
    
    [Ainvt,posdef]          = chol(toeplitz(dRESac_adj)); 
    p1                      = size(Ainvt,1); % this is basically posdef - 1, but we'll keep it as Keith Worsely's code. 
    A                       = inv(Ainvt'); 
    
    sqrtmVhalf              = toeplitz([A(p1,p1:-1:1) zeros(1,T-p1)],zeros(1,T)); 
    sqrtmVhalf(1:p1,1:p1)   = A;
    
%     Coradj_pix=squeeze(rho_vol(pix,slice,:));
%     [Ainvt posdef]=chol(toeplitz([1 Coradj_pix']));
%     tic;
%     B=ones(T-p1,1)*A(p1,:);
%     V00=spdiags(B,1:p1,T-p1,T); 
%     toc
end