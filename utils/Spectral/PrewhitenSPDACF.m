function ts0 = PrewhitenSPDACF(Yt,ARorder,Fs)

V    = size(Yt,2); 
T    = size(Yt,1); %this can be replaced with S = sum(win.^2) if/when the windowing is used
nfft = T; % keep it like this, if the one sided hacked, this will be handy


% [ARyEst,ARv] = AR_YW(Yt,ARorder); % sanity check
% h            = freqz(1,[1 -ARyEst'],nfft,'whole',Fs); 

% [ARyEst,ARv] = aryule(Yt,ARorder); % sanity check
% 
% h = zeros(T,V);
% for vv = 1:V
%     h(:,vv) = freqz(1,ARyEst(vv,:),nfft,'whole',Fs);
% end

% Yt      = Yt-mean(Yt);
% Yt      = Yt./std(Yt); 

acf          = AC_fft(Yt,T);
acf          = acf'; 
acf          = [1;-acf(2:ARorder)];
acf(1:5)
h = freqz(1,acf,nfft,'whole',Fs);

%acf(2:end,:) = -acf(2:end,:);
%h = freqz(1,acf,895,'whole',1);
  ARv          = 1; % Or the Bartlett's var(acf) = (1+2.*sum(acf(1:round(T/4),:).^2))./T ?
%  acf          = acf(1:ARorder,:); 
%  acf          = [acf; zeros(T-ARorder,V)];
% % size(acf)
%  acf_fft = fft(acf,nfft,1);
%  h       = 1./acf_fft(1:nfft,:); 

%%% SANITY CHECK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % freqz gives us inverse of the fft of the ARYW parametes
% ARorder = 10;
% [ARyEst,ARv] = aryule(Y(:,5),ARorder); % sanity check
% g = freqz(1,ARyEst,895,'whole',1);
% plot(real(g))
% hold on; 
% gg = fft([ARyEst zeros(1,895-ARorder)],895);
% plot(real(1./gg))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
figure; hold on; 
plot(rawpxx)
plot(Pxx0)

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