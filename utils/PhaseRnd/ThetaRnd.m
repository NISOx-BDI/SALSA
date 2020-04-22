function sY = ThetaRnd(Y,T)
%
% Generate time series using phase randomisation. 
% 
% Y: should be a TxV matrix  
% T: Number of datapoints. Must be odd
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
if ismember(size(Y),T); erorr('ThetaRnd:: Check input.'); end;
if rem(T,2)==0; erorr('ThetaRnd:: time series should be odd.'); end;
if size(Y,1)~=T; Y=Y'; end; 

[T,nts] = size(Y);

% Get parameters
len_ser     = (T-1)/2;
RightSide   = 2:len_ser+1; 
LeftSide    = len_ser+2:T;

% Fourier transform of the original dataset
fft_recblk = fft(Y);

% Generate random phases between 0-2\pi
RndPhase         = exp(2*pi*1i*rand(len_ser,1));
RndPhase_RIGHT   = repmat(RndPhase,1,nts);
% Phase should be symmetrised, so flip and use conjugate on the complex numbers
RndPhase_LEFT    = conj(flipud(RndPhase_RIGHT));

% Randomize all the time series simultaneously
fft_recblk_surr              = fft_recblk;
fft_recblk_surr(RightSide,:) = fft_recblk(RightSide,:).*RndPhase_RIGHT;
fft_recblk_surr(LeftSide,:)  = fft_recblk(LeftSide,:).*RndPhase_LEFT;

% Inverse transform
sY= real(ifft(fft_recblk_surr));

