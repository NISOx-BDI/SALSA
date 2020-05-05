function ts0 = PrewhitenSPD(Yt,ARorder,Fs)

V    = size(Yt,2); 
T    = size(Yt,1); %this can be replaced with S = sum(win.^2) if/when the windowing is used
nfft = T; % keep it like this, if the one sided hacked, this will be handy


% [ARyEst,ARv] = AR_YW(Yt,ARorder); % sanity check
% h            = freqz(1,[1 -ARyEst'],nfft,'whole',Fs); 

[ARyEst,ARv] = aryule(Yt,ARorder); % sanity check

h = zeros(T,V);
for vv = 1:V
    h(:,vv) = freqz(1,ARyEst(vv,:),nfft,'whole',Fs);
end

Sxx     = bsxfun(@times,ARv,abs(h).^2);   % make it power & multiply by variance of AR 
Pxx0    = Sxx./Fs; % this is theoritical periodgram of the AR/ARMA process
%%%% two sided get the periodogram 
xdft    = fft(Yt,nfft,1);
psdx    = abs(xdft).^2; %xdft.*conj(xdft); % psdx = abs(xdft).^2
SF      = (1/(Fs*T)); % to get a two-sided, the power should be doubled. 
rawpxx  = SF * psdx;

% these are only useful for plotting & sanity check
%freq   = 0:Fs/nfft:Fs; 
%freq   = freq(1:end-1); 

Wspec   = rawpxx./Pxx0; % flatten the periodogram

Sxx1    = Wspec; % remove the factor
Sxx1    = Sxx1.*T.*Fs;
ff      = sqrt(Sxx1);

phase0 = zeros(T,V);
for vv = 1:V
    phase0(:,vv) = phase(xdft(:,vv));
end

ts0    = real(ifft(ff.*exp(1j*phase0),nfft));

% The windowing is shouldn't be done here. 
%win  = rectwin(T); 
%Yt   = Yt.*win; % Apply windowing here (periodogram does it on it's own)