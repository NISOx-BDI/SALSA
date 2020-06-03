function [CPS,stat,CPZ] = CPSUnivar(Resid,X)
% Resid : resiuldas [Time x Voxel]
% X     : design [Time x P]
% 
% TEN & SA, Ox, 2020

nScan = size(Resid,1);
nVar = rank(X);
nVox = size(Resid,2); 

Spec = abs(fft(Resid,nScan-nVar)).^2;

% Select first half  
if rem(nScan-nVar,2)         % nfft odd
    select   = (2:(nScan-nVar+1)/2)';  %-don't include DC or Nyquist components 
else
    select   = (2:(nScan-nVar)/2)';
end

%freq  = (select-1)*2/(nScan-nVar);
power = [2*Spec(select,:)];


CPseries = cumsum(power);
Flength = size(CPseries,1);
CPS	= CPseries./(repmat(sum(power),Flength,1));

%-setup the upper and lower statistic for the Cumulative periodogram
MaxA1 = abs(CPS(1:Flength-1,:)-...
    (repmat([1:Flength-1]',1,nVox))/(Flength));
MaxB1 = abs(CPS(1:Flength-1,:)-...
    (repmat([1:Flength-1]',1,nVox)-1)/(Flength));

%-Calculate the Cumulative periodogram test statistic
stat  = max(max(MaxA1,MaxB1));
CPZ  = stat*(sqrt(Flength-1)+0.12+0.11/sqrt(Flength-1));
% p    = 1-sqrt(2*pi)*(sum(exp(-((2*repmat([1:30]',1,nVox)-1)*pi).^2./(8*repmat(CPZ,size([1:30],2),1).^2))))./CPZ;
% p    = p+(p==0)*eps;
% p    = -log10(p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% OLD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [CPS,stat,CPZ] = CPSUnivar(Resid,X)
% 
% nScan = numel(Resid);
% nVar = rank(X);
% nVox = 1; 
% 
% Spec = abs(fft(Resid,nScan-nVar)).^2;
% 
% % Select first half  
% if rem(nScan-nVar,2)         % nfft odd
%     select   = (2:(nScan-nVar+1)/2)';  %-don't include DC or Nyquist components 
% else
%     select   = (2:(nScan-nVar)/2)';
% end
% 
% freq  = (select-1)*2/(nScan-nVar);
% power = [2*Spec(select)];
% 
% CPseries = cumsum(power);
% Flength = size(CPseries,1);
% CPS	= CPseries./(repmat(sum(power),Flength,1));
% 
% %-setup the upper and lower statistic for the Cumulative periodogram
% MaxA1 = abs(CPS(1:Flength-1,:)-...
%     (repmat([1:Flength-1]',1,nVox))/(Flength));
% MaxB1 = abs(CPS(1:Flength-1,:)-...
%     (repmat([1:Flength-1]',1,nVox)-1)/(Flength));
% 
% %-Calculate the Cumulative periodogram test statistic
% stat  = max(max(MaxA1,MaxB1));
% CPZ  = stat*(sqrt(Flength-1)+0.12+0.11/sqrt(Flength-1));
% % p    = 1-sqrt(2*pi)*(sum(exp(-((2*repmat([1:30]',1,nVox)-1)*pi).^2./(8*repmat(CPZ,size([1:30],2),1).^2))))./CPZ;
% % p    = p+(p==0)*eps;
% % p    = -log10(p);