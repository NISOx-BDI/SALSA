clear

T    = 500;
TR   = 3; 
Fs   = 1/TR;
nRlz = 1; 
%simulate data using filter
a    = [1 -0.9];%ARcoeffs
b    = 1 ;%MAcoeffs

for i = 1:nRlz
    e = randn(T,1); % generate gaussian white noise
    Yt(:,i) = filter(b,a,e); % generate y
end

%Yt = Yt./std(Yt);
Yt = Yt - mean(Yt);

% --- ESTIMATE DENSITY -- 
nfft=T;
[ARyEst,ARv] = aryule(Yt,6); % sanity check
[h,w]        = freqz(1,ARyEst,nfft,'whole',Fs);
Sxx          = bsxfun(@times,ARv,abs(h).^2);   % make it power multiply by variance of AR 

%one-sided
%select       = 1:nfft/2+1;    % EVEN
%Sxx_unscaled = Sxx(select); % Take only [0,pi] or [0,pi)
%Sxx          = [Sxx_unscaled(1); 2*Sxx_unscaled(2:end-1); Sxx_unscaled(end)];
%w             = w(select);
Pxx0          = Sxx./Fs;
%-----

[rawpxx,wp] = periodogram(Yt,rectwin(T),nfft,Fs,'two-sided');

[rawpxx1s,wp1s] = periodogram(Yt,rectwin(T),nfft,Fs);

%%%% two sided get the periodogram 
win = rectwin(T); 
xw   = Yt.*win; % Apply windowing here (periodogram does it on it's own)
xdft = fft(Yt,nfft);
psdx = abs(xdft).^2; %xdft.*conj(xdft); % psdx = abs(xdft).^2
S    = sum(win.^2);
SF   = (1/(Fs*S));
rawpxx0 = SF * psdx;
freq = 0:Fs/nfft:Fs; freq = freq(1:end-1); 



Wspec   = rawpxx./Pxx0;
%Wspec    = rawpxx;

Sxx1 = Wspec; % remove the factor
% hfreq = Sxx1(end); % middle of the whole spectrum
% lfreq = Sxx1(1);   % the first & the last point of the spect.
% Sxx1 = [lfreq; Sxx1(2:end-1)/2; hfreq]; % make the half of one sided
% Sxx1 = [lfreq;Sxx1(2:end-1); flip(Sxx1(2:end-1));lfreq]; % make the spec double sided
Sxx1 = Sxx1*T*Fs;
ff   = sqrt(Sxx1);
 
xdft   = fft(Yt,nfft);
phase0 = phase(xdft);
ts0    = real(ifft(ff.*exp(1j*phase0),nfft));

figure; hold on; plot(ts0); plot(Yt)

method = @aryule;
[Pxx,f] = arspectra00(method,Yt,1,nfft,Fs);

%[Pxx,w,units] = computepsd(Sxx,w); 

% draw the spectrum of the simulated time series
%[xx,yy] = DrawMeSpectrum(Yt,TR,-1);

%figure; 
%hold on 
%plot(xx,mean(yy,2),'color','r','linewidth',1.3)

% parametric
[sdc,wc] = pcov(Yt,1,nfft,Fs,'two-sided'); % estimate spectral density by fitting AR(8)
[sdy,wy] = pyulear(Yt,1,nfft,Fs,'two-sided'); % estimate spectral density by fitting AR(8) using the Yule?walker equations

figure; 
hold on; grid on; 
plot(wp,rawpxx,wc,sdc,wy,sdy,wy,Pxx0,wy,Wspec);
legend('raw periodogram','COV AR','Yule-Walker','My YW','PW')%,'population density');

% % nonparametric
%[sdp,wp] = periodogram(Yt,[],T,'onesided'); % estimate using unsmoothed
% % periodogram
% rT = round(sqrt(T))*3;
% [sdw,ww] = pwelch(Yt,parzenwin(rT),rT-1,T,'onesided'); % smoothed periodogram

%bartlett weighted covariances
% [c,lags]=xcov(Yt,'biased');
% t = -(T-1):(T-1);
% weight = 1-abs(t')/rT;
% weight(abs(t')>rT) = 0;
% wb=ww;
% for j=1:length(wb)
%     sdb(j) = sum(c.*weight.*exp(-i*wb(j)*(-(T-1):(T-1))'));
% end
% sdb = sdb / sqrt(2*pi);                                                                                     

% %population density
% w=ww;
% h = freqz(b,a,w);
% sd = abs(h).^2./sqrt(2*pi);

%figure;
%plot(wp,sdp,wc,sdc,ww,sdw,wb,sdb,wy,sdy,w,sd);
%legend('raw periodogram','parametric AR','smoothed periodogram', ...
%'bartlett weighted cov','Yule?Walker','population density');

%plot(wp,sdp,wc,sdc,ww,sdw,wy,sdy)%,w,sd)%,wb,sdy);

%legend('raw periodogram','parametric AR','smoothed periodogram','Yule-Walker')%,'population density');


% Fs=1;
% xdft = fft(Yt);
% xdft = xdft(1:T/2+1);
% psdx = (1/(Fs*T)) * abs(xdft).^2;
% psdx(2:end-1) = 2*psdx(2:end-1);
% freq = 0:Fs/length(Yt):Fs/2;
% figure; hold on; 
% plot(freq,psdx)
% plot(wp,sdp)
% grid on
% title('Periodogram Using FFT')
% xlabel('Frequency (Hz)')
% ylabel('Power/Frequency (dB/Hz)')

% [xx,f] = pyulear(Yt,1,1/TR);
% 
% plot(f,xx)