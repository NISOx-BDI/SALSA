function ts0 = Prewhiten2x2SPD(Yt2,Yt1,ARorder,Fs)

V1    = size(Yt1,2); 
T1    = size(Yt1,1); 


V2    = size(Yt2,2); 
T2    = size(Yt2,1); %this can be replaced with S = sum(win.^2) if/when the windowing is used
nfft = T2; % keep it like this, if the one sided hacked, this will be handy



% [ARyEst,ARv] = AR_YW(Yt,ARorder); % sanity check
% h            = freqz(1,[1 -ARyEst'],nfft,'whole',Fs); 

[ARyEst,ARv] = aryule(Yt1,ARorder); % sanity check

h = zeros(T1,V1);
for vv = 1:V1
    h(:,vv) = freqz(1,ARyEst(vv,:),nfft,'whole',Fs);
end

Sxx     = bsxfun(@times,ARv,abs(h).^2);   % make it power & multiply by variance of AR 
Pxx0    = Sxx./Fs; % this is theoritical periodgram of the AR/ARMA process
%%%% two sided get the periodogram 
xdft    = fft(Yt2,nfft,1);
psdx    = abs(xdft).^2; %xdft.*conj(xdft); % psdx = abs(xdft).^2
SF      = (1/(Fs*T2)); % to get a two-sided, the power should be doubled. 
rawpxx  = SF * psdx;

% these are only useful for plotting & sanity check
%freq   = 0:Fs/nfft:Fs; 
%freq   = freq(1:end-1); 

Wspec   = rawpxx./Pxx0; % flatten the periodogram

Sxx1    = Wspec; % remove the factor
Sxx1    = Sxx1.*T2.*Fs;
ff      = sqrt(Sxx1);

phase0 = zeros(T2,V2);
for vv = 1:V2
    phase0(:,vv) = phase(xdft(:,vv));
end

ts0    = real(ifft(ff.*exp(1j*phase0),nfft));






% function ts0 = Prewhiten2SPD(ED,Yt,ARorder,Fs)
% 
% if isempty(ARorder); ARorder = 4; end;  
% 
% T    = numel(Yt); %this can be replaced with S = sum(win.^2) if/when the windowing is used
% nfft = T; % keep it like this, if the one sided hacked, this will be handy
% 
% %[ARyEst,ARv] = aryule(Yt,6); % sanity check, note the negative AR par
% [ARyEst,ARv] = AR_YW(Yt,ARorder);
% h            = freqz(1,[1 -ARyEst'],nfft,'whole',Fs); 
% Sxx          = bsxfun(@times,ARv,abs(h).^2);   % make it power & multiply by variance of AR 
% Pxx0         = Sxx./Fs; % this is theoritical periodgram of the AR/ARMA process
% 
% %%%% two sided get the periodogram 
% xdft   = fft(ED,nfft);
% psdx   = abs(xdft).^2; %xdft.*conj(xdft); % psdx = abs(xdft).^2
% SF     = (1/(Fs*T)); % to get a two-sided, the power should be doubled. 
% rawpxx = SF * psdx;
% 
% % these are only useful for plotting & sanity check
% %freq   = 0:Fs/nfft:Fs; 
% %freq   = freq(1:end-1); 
% 
% Wspec  = rawpxx./Pxx0; % flatten the periodogram
% 
% Sxx1   = Wspec; % remove the factor
% Sxx1   = Sxx1.*T.*Fs;
% ff     = sqrt(Sxx1);
%  
% phase0 = phase(xdft);
% ts0    = real(ifft(ff.*exp(1j.*phase0),nfft));
% 
% % The windowing is shouldn't be done here. 
% %win  = rectwin(T); 
% %Yt   = Yt.*win; % Apply windowing here (periodogram does it on it's own)