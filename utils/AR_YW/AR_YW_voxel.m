
function rho = AR_YW_voxel(Y,T,p)
% rho = AR_YW_voxel(Y,p)
%
%%% INPUTS
% Y : time series -- it can only be a vector/column for time being
% T : length of the time series -- for sanity check -- intrger
% p : order -- integer
%
%%% OUTPUTS
% rho : AR parameters
%
% SA, Ox, 2020

    ac  = AC_fft(Y,T);    
    %ac(round(2*sqrt(T)):end) = 0; % this might mess up stuff!
    [nv,nl] = size(ac);
    rho     = zeros(nv,p);
    
    for vi = 1:nv
        if ~mod(vi,10000); disp(['AR_YW_voxel:: on voxel ' num2str(vi)]); end; 
        R         = toeplitz(ac(vi,1:p));
        r         = ac(vi,2:p+1);    
        rho(vi,:) = pinv(R)*r';
    end
end

function [xAC,CI,ACOV]=AC_fft(Y,L,varargin)
%[xAC]=AC_fft(Y,T,varargin)
% Super fast full-lag AC calculation of multi-dimention matrices. The
% function exploits fft to estimate the autocorrelations. 
%
%%%%INPUTS
%   Y:      A matrix of size IxT comprised of I time series of T length.
%   L:      Time series length
%   
%   To get a double sided AC, add 'two-sided' as an input.
%
%   NB! The function includes the 0lag zero (i.e. 1) to the output. 
%%%%OUTPUTS
%   xAC:    IxT-1 matrix of full-lag autocorrelations. If 'two-sided'
%           switched, then the matrix is Ix2(T-1).
%   CI :    95% Confidence Intervals of AFC.
%
%_________________________________________________________________________
% Soroosh Afyouni, University of Oxford, 2017
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________

if size(Y,2)~=L
    Y=Y';
end

if size(Y,2)~=L
    error('Use IxT or TxI input form.')
end

Y=Y-mean(Y,2); 
%only works on >2016 Matlabs, but faster!
% if <2016, use Y=Y-repmat(mean(Y,2),1,L) instead.

nfft    = 2.^nextpow2(2*L-1); %zero-pad the hell out!
yfft    = fft(Y,nfft,2); %be careful with the dimensions

ACOV = real(ifft(yfft.*conj(yfft),[],2));

ACOV = ACOV(:,1:L);

%xAC = ACOV;
xAC  = ACOV./sum(abs(Y).^2,2); %normalise the COVs

if sum(strcmpi(varargin,'two-sided')) %two sided is just mirrored, AC func is symmetric
   xAC  = [xAC(:,end-L+2:end) ; xAC(:,1:L)];
else
    xAC  = xAC(:,1:L);
end

bnd=(sqrt(2)*erfinv(0.95))./sqrt(L); %assumes normality for AC
CI=[-bnd bnd];

end