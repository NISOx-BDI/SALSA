function [sqrtmVhalf,spdflag] = ACF_ResPWm(autocov,AROrd,BiasAdj,tflag)
%
% dRESacov : autocovariance 
% Mord : Order required 
% BiasAdj: bias adjustment matrix 
%
% SA, Ox, 2020
%

    autocov = autocov(:); %make everything a column!
    if nargin < 3 || isempty(BiasAdj); BiasAdj = eye(AROrd+1); end
    if nargin < 2; AROrd = 2*round(sqrt(T)); end
    
    autocov              =  autocov(:);    
    T                    = numel(autocov); % full ACF has one element less than the time series    
    dRESac_adj           = (BiasAdj*autocov(1:AROrd+1));
    dRESac_adj           = dRESac_adj./dRESac_adj(1); % make a auto*correlation*
    
    if tflag
       %disp('Tukey taper.')
       acfdRES_tt       = [1 TukeyTaperMe(dRESac_adj(2:end),T-1,AROrd)];
    else
       %disp('NO tukey taper')
       acfdRES_tt        = [dRESac_adj;zeros(T-AROrd-1,1)];
    end
    
    
    
    %-------- sanity check ---------
%     acfdRES_tt(1:AROrd+5)
%     size(acfdRES_tt)
%     figure; hold on; grid on; box on; 
%     plot(dRESac_adj(1:AROrd))
%     plot(acfdRES_tt(1:AROrd+5))
%     xlim([0 50]);
    %-------------------------------
    
    
    ACMat                = toeplitz(acfdRES_tt);  
    [sqrtmVhalf,spdflag] = CholWhiten(ACMat); 
end

function tt_ts = TukeyTaperMe(acs,T,M)
% acs MUST be a column. 
%SA, Ox, 2018
    M          = round(M);
    tt_ts      = zeros(1,T);
    tt_ts(1:M) = (1+cos([1:M]'.*pi./M))./2.*acs(1:M);    
end