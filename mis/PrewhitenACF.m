function pwY = PrewhitenACF(Y,AROrd,BiasAdj,tflag)
% pwY = PrewhitenACF(Y,AROrd,BiasAdj,tflag)
% Prewhiten a given time series using autocorrelation function. 
% If autocorrelation matrix is not PSD, it finds the nearest PSD. 
%
%%%%% INPUTS:
% Y      : The time series. Must be TxV where T is length of time series
%          and V is number of ROIs/voxels.
% Mord   : Order of ACF , if not specified, it will be set to sqrt(T) (Woolrich et al 2001)
% BiasAdj: bias adjustment matrix. If not specified it will be eye(AROrd+1)
% tflag  : If switched, then the ACF is regularised by a Tukey taper.
%          Defualt is 1
% 
%%%% OUTPUT:
% 
% pwY : prewhitened time series
% 
%%%% EXAMPLE:
% Y = filter(1,[1 -0.6],randn(500,10)); % generate 10 time series of AR(1)=0.6
% autocorr(Y(:,2)) % check the ACF
% pwY = PrewhitenACF(Y); % prewhiten Y
% autocorr(pwY(:,2)) % check the ACF of prewhitened time series [should ideally be flat]
%
%_________________________________________________________________________
% Soroosh Afyouni, Ox, 2020
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________

    [T,V]           = size(Y);

    if nargin < 2 || isempty(AROrd);   AROrd = round(sqrt(T)); end 
    if nargin < 3 || isempty(BiasAdj); BiasAdj = eye(AROrd+1); end
    if nargin < 4 || isempty(tflag);   tflag = 1; end; 
    
    [~,~,ACOV]   = AC_fft(Y,T); %faster if ACF is calculated all at once here. 
    for iv = 1:V
        autocov = ACOV(iv,:)';
        %autocov = autocov(:); %make everything a column!

        autocov              =  autocov(:);    
        dRESac_adj           = (BiasAdj*autocov(1:AROrd+1));
        dRESac_adj           = dRESac_adj./dRESac_adj(1); % make an auto*correlation*

        if tflag
           %disp('Tukey taper.')
           acf_tt       = [1 TukeyTaperMe(dRESac_adj(2:end),T-1,AROrd)];
        else
           %disp('NO tukey taper')
           acf_tt        = [dRESac_adj;zeros(T-AROrd-1,1)];
        end

        ACMat                = toeplitz(acf_tt);  
        [Vhalf,spdflag] = CholWhiten(ACMat); 
        pwY = Vhalf*Y; 
        
    %-------- sanity check ---------
    %     acfdRES_tt(1:AROrd+5)
    %     size(acfdRES_tt)
    %     figure; hold on; grid on; box on; 
    %     plot(dRESac_adj(1:AROrd))
    %     plot(acfdRES_tt(1:AROrd+5))
    %     xlim([0 50]);
    %-------------------------------
            
    end
    
end

function tt_ts = TukeyTaperMe(acs,T,M)
% acs MUST be a column. 
%SA, Ox, 2018
    M          = round(M);
    tt_ts      = zeros(1,T);
    tt_ts(1:M) = (1+cos([1:M]'.*pi./M))./2.*acs(1:M);    
end

function [Vhalf,spdflag] = CholWhiten(COVmat)
% [W,spdflag] = CholWhiten(COV)
% Whiten with Cholesky decomposition & inverse
% SA, OX, 2020
%
  [R,spdflag]   = chol(COVmat);
  if spdflag
      disp('PrewhitenACF:: CholWhiten:: used nearestSPD.')
      spdCOVmat    = nearestSPD(COVmat);
      R            = chol(spdCOVmat);
  end
  Vhalf         = inv(R');
end

function Ahat = nearestSPD(A)
% From Higham: "The nearest symmetric positive semidefinite matrix in the
% Frobenius norm to an arbitrary real matrix A is shown to be (B + H)/2,
% where H is the symmetric polar factor of B=(A + A')/2."
%
% http://www.sciencedirect.com/science/article/pii/0024379588902236
%-------------------------------------------------------------------------
%Copyright (c) 2013, John D'Errico 
% All rights reserved.
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% * Redistributions of source code must retain the above copyright 
% notice, this list of conditions and the following disclaimer. 
% * Redistributions in binary form must reproduce the above copyright 
% notice, this list of conditions and the following disclaimer in 
% the documentation and/or other materials provided with the distribution

if nargin ~= 1
  error('Exactly one argument must be provided.')
end
[r,c] = size(A);
if r ~= c
  error('A must be a square matrix.')
elseif (r == 1) && (A <= 0)
  Ahat = eps;
  return
end
B = (A + A')/2;
[~,Sigma,V] = svd(B);
H    = V*Sigma*V';
Ahat = (B+H)/2;
Ahat = (Ahat + Ahat')/2;
p    = 1;
k    = 0;
while p ~= 0
  [~,p]  = chol(Ahat);
  %disp('nearestSPD:: Still not there!')
  k = k + 1;
  if p  ~= 0
    mineig = min(eig(Ahat));
    Ahat   = Ahat + (-mineig*k.^2 + eps(mineig))*eye(size(A));
  end
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