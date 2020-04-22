
% % ------- A simple cosine wave
% 
% fc      = 10;%frequency of the cosine wave
% fs      = 32*fc;%sampling frequency with oversampling factor 32
% t       = 0:1/fs:2-1/fs;    %2 seconds duration
% phi     = 0; % no phase shift for time being
% xorig   = cos(2*pi*fc*t+phi);
% xorig   = xorig./std(xorig);
% 
% figure; 
% hold on;
% plot(xorig(1:end-1)); %plot the signal
% 
% XX = phaseran(xorig',10);
% for i = 1:10
%     plot(XX(:,1,i)); 
% end
% 
% figure; 
% hold on; 
% [xs,ys] = DrawMeSpectrum(xorig,1/fs);
% plot(xs,ys)
% [xxs,yys] = DrawMeSpectrum(squeeze(XX(:,1,:)),1/fs);
% plot(xxs,yys)
% 
% % ------- 1/f noise
% % copied from: https://uk.mathworks.com/matlabcentral/answers/344786-generation-of-1-f-noise-using-matlab
% fv = linspace(0, 1, 20);                                % Normalised Frequencies
% a = 1./(1 + fv*2);                                      % Amplitudes Of ?1/f?
% b = firls(42, fv, a);                                   % Filter Numerator Coefficients                                     % Filter Bode Plot
% ns = rand(1, T);
% xorig = filtfilt(b, 1, ns);                             % Create ?1/f? Noise
% xorig = xorig./std(xorig);
% 
% figure; 
% hold on;
% plot(xorig(1:end-1)); %plot the signal
% 
% XX = phaseran(xorig',10);
% for i = 1:10
%     plot(XX(:,1,i)); 
% end
% 
% figure; 
% hold on; 
% [xs,ys] = DrawMeSpectrum(xorig,1/fs);
% plot(xs,ys)
% [xxs,yys] = DrawMeSpectrum(squeeze(XX(:,1,:)),1/fs);
% plot(xxs,yys)

% -------- Real Data ----
PATH2AUX='/Users/sorooshafyouni/Home/GitClone/FILM2';
addpath([PATH2AUX '/utils/Trend'])
addpath([PATH2AUX '/mis'])

SubID     = 'A00027167';
SesID     = 'DS2'; 
TR        = 0.645; 
disp('=======================================')
PATH2AUX='/Users/sorooshafyouni/Home/GitClone/FILM2';
addpath([PATH2AUX '/utils/Trend'])
addpath('/Users/sorooshafyouni/Home/matlab/spm12')
%12487
disp('=====SET UP PATHS =============================')
%Raw Images (MMP feat output)
Path2ImgRaw=[PATH2AUX '/ExampleData/R.mpp'];
%Path2ImgDir = [Path2ImgRaw '/sub-' SubID '/ses-' SesID '/sub-' SubID '_ses-' SesID '_task-rest_acq-' num2str(TR*1000) '_bold.mpp'];
Path2ImgDir = ['/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/R_test/sub-' SubID '_ses-' SesID '_task-rest_acq-645_bold_mpp'];
Path2Img    = [Path2ImgDir '/prefiltered_func_data_bet.nii'];
xorig = CleanNIFTI_spm(Path2Img,'demean');

xorig = xorig(:,1:end-1); % odd numbers
xorig = xorig'; % TxV

% generat surrogate data
XX = ThetaRnd(xorig,size(xorig,1));

%-- Sanity check to ensure the acfs are pretty similar
% xorig_ACL   = sum(AC_fft(xorig,size(xorig,1)).^2);
% XX_ACL      = sum(AC_fft(XX,size(xorig,1)).^2);
% scatter(XX_ACL,xorig_ACL,'.')

%-- Sanity check to ensure the spectrums are pretty similar

figure; 
subplot(4,1,1)

hold on; 
[xs,ys] = DrawMeSpectrum(xorig,TR);
plot(xs,mean(ys,2),'LineWidth',2,'color',[.5 .5 .5])
ylabel('Power'); xlabel('Freq')

[xxs,yys] = DrawMeSpectrum(XX,TR);
plot(xxs,mean(yys,2),'LineWidth',2,'color',[1 0 0])
ylabel('Power'); xlabel('Freq')
line([0 max(xxs)],[1 1],'color','k','linestyle','-.')

% -- HPF filter to detrend
hpf_filter = hp_fsl(size(xorig,1),100,TR);

% filter the original data
dX = hpf_filter*xorig;
% filter the surrogate data
dXX = hpf_filter*XX;

subplot(4,1,2)
hold on; 
[dxs,dys] = DrawMeSpectrum(dX,TR);
plot(dxs,mean(dys,2),'LineWidth',2,'color',[.5 .5 .5])
ylabel('Power'); xlabel('Freq')

[dxxs,dyys] = DrawMeSpectrum(dXX,TR);
plot(dxxs,mean(dyys,2),'LineWidth',2,'color',[1 0 0])
ylabel('Power'); xlabel('Freq')
line([0 max(xxs)],[1 1],'color','k','linestyle','-.')

% --- polynomials to detrend

pnum = 1+floor(TR.*size(xorig,1)/150);
dX = multpolyfit(repmat(1:size(xorig,1),size(xorig,2),1)',xorig,size(xorig,1),pnum);
dXX = multpolyfit(repmat(1:size(xorig,1),size(xorig,2),1)',XX,size(xorig,1),pnum);

subplot(4,1,3)
hold on; 
[dxs,dys] = DrawMeSpectrum(dX',TR);
plot(dxs,mean(dys,2),'LineWidth',2,'color',[.5 .5 .5])
ylabel('Power'); xlabel('Freq')

[dxxs,dyys] = DrawMeSpectrum(dXX',TR);
plot(dxxs,mean(dyys,2),'LineWidth',2,'color',[1 0 0])
ylabel('Power'); xlabel('Freq')
line([0 max(xxs)],[1 1],'color','k','linestyle','-.')

% N=256; %FFT size
% X = 1/N*fftshift(fft(x,N));%N-point complex DFT
% 
% 
% df=fs/N; %frequency resolution
% sampleIndex = -N/2:N/2-1; %ordered index for FFT plot
% f=sampleIndex*df; %x-axis index converted to ordered frequencies
% stem(f,abs(X)); %magnitudes vs frequencies
% xlabel('f (Hz)'); ylabel('|X(k)|');
% 
% phase1 = unwrap(angle(X).*256);
% 
% phase=atan2(imag(X),real(X))*180/pi; %phase information
% plot(f,phase); %phase vs frequencies
% 
% X(1:5)
% atan2(imag(X(1:5)),real(X(1:5)))
% 
% X2=X;%store the FFT results in another array
% %detect noise (very small numbers (eps)) and ignore them
% threshold = max(abs(X))/10000; %tolerance threshold
% X2(abs(X)<threshold) = 0; %maskout values that are below the threshold
% 
% phase=atan2(imag(X2),real(X2))*180/pi; %phase information
% plot(f,phase); %phase vs frequencies