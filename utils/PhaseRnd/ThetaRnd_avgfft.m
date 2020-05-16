function sY = ThetaRnd_avgfft(Y_dft,T,seedid)
%
% Generate surrogate time series using phase randomisation. 
% 
% Y        : should be a TxV matrix  
% T        : Number of datapoints. Must be odd
% corrflag : if unswitched, then it break the cross-corr structures [default: 1]
%
% Reference:
% 
% https://uk.mathworks.com/matlabcentral/fileexchange/32621-phase-randomization
% 
% Theiler, James, et al. Testing for nonlinearity in time series: the 
% method of surrogate data. No. LA-UR-91-3343; CONF-9108181-1. Los Alamos 
% National Lab., NM (United States), 1991.
% 
% Prichard, D., Theiler, J. Generating Surrogate Data for Time Series
% with Several Simultaneously Measured Variables (1994)
% Physical Review Letters, Vol 73, Number 7
% Get parameters
%
% SA, Ox, 2020
% 


% check the input 
%if ~any(ismember(size(Y),T)); error('ThetaRnd:: Check input.'); end;
if nargin==3; rand('seed',seedid);     end;


%[T,nts]     = size(Y);
nts = 1;

if ~mod(T,2) % if the time series are even, exclude the end point from analysis
    error('The spectrum should be an odd vector in length')
end

% Get parameters
len_ser     = (T-1)/2;
RHS_idx     = 2:len_ser+1; % Is this because the first phase is zero?!
LHS_idx     = len_ser+2:T;

% DFT
%Y_dft     = fft(Y);

% Generate random phases between 0-2\pi

RndPhase         = exp(2*pi*1i*rand([len_ser,1]));
RndPhase_RIGHT   = repmat(RndPhase,1,nts);
    
% Phase should be symmetrised, so flip and use conjugate on the complex numbers
RndPhase_LEFT    = conj(flipud(RndPhase_RIGHT));

% Use random phases. If needed to break the cross-corr structure, use
% random order of phases, probablu usin permute.m?
Yr_dft             = Y_dft;
Yr_dft(RHS_idx,:)  = Y_dft(RHS_idx,:).*RndPhase_RIGHT;
Yr_dft(LHS_idx,:)  = Y_dft(LHS_idx,:).*RndPhase_LEFT;

% ifft
sY= [real(ifft(Yr_dft))];

