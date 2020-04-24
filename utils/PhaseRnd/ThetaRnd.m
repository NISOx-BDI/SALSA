function sY = ThetaRnd(Y,T,corrflag)
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
if ismember(size(Y),T); erorr('ThetaRnd:: Check input.');       end;
if rem(T,2)==0; erorr('ThetaRnd:: time series should be odd.'); end;
if size(Y,1)~=T; Y=Y';                                          end; 
if nargin<3; corrflag = 1;                                      end; 

[T,nts]     = size(Y);

% Get parameters
len_ser     = (T-1)/2;
RHS_idx     = 2:len_ser+1; % Is this because the first phase is zero?!
LHS_idx     = len_ser+2:T;

% DFT
Y_dft     = fft(Y);

% Generate random phases between 0-2\pi
if corrflag
    RndPhase         = exp(2*pi*1i*rand([len_ser,1]));
    RndPhase_RIGHT   = repmat(RndPhase,1,nts);
else
    RndPhase_RIGHT   = exp(2.*pi.*1i.*rand([len_ser,nts]));
end
    
% Phase should be symmetrised, so flip and use conjugate on the complex numbers
RndPhase_LEFT    = conj(flipud(RndPhase_RIGHT));

% Use random phases. If needed to break the cross-corr structure, use
% random order of phases, probablu usin permute.m?
Yr_dft             = Y_dft;
Yr_dft(RHS_idx,:)  = Y_dft(RHS_idx,:).*RndPhase_RIGHT;
Yr_dft(LHS_idx,:)  = Y_dft(LHS_idx,:).*RndPhase_LEFT;

% ifft
sY= real(ifft(Yr_dft));

