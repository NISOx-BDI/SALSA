clear

addpath(genpath('/Users/sorooshafyouni/Home/GitClone/xDF'))
addpath('/Users/sorooshafyouni/Home/matlab/spm12')
addpath('/Users/sorooshafyouni/Home/BCF/ARMA/AR_YW')

TR = 2; 

load(['/Users/sorooshafyouni/Home/PREW/SimulatefMRI/S_1D/ts_TR' num2str(TR*1000) '_RHOt48_LF128.mat'])
[Nt,Nl] = size(ts);

%%% EXPERIMENTAL DESIGN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stick design
% Nsticks = 30;
% onsets=randperm(Nt); 
% DD=zeros(1,Nt); 
% DD(onsets(1:Nsticks))=1;

% boxcar design
bcduration=20;
BCWidth = bcduration/TR;
oneSTIM = [zeros(1, BCWidth), ones(1, BCWidth)];
numberSTIM = Nl/numel(oneSTIM);
DD = repmat(oneSTIM, [1, numberSTIM]);

hrf=spm_hrf(TR); 
conv_dsgn=conv(DD,hrf); 
Xraw=conv_dsgn(1:Nl)';

X = Xraw - mean(Xraw); 
design_eff=1/trace(inv(X'*X));

%%% INDUCE TASK? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for it = 1:Nt
%    b_hat(it)=X\ts(it,:)';
%end

%%% DETEREND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ACts = AC_fft(ts,Nl);
HPFfilter=hp_fsl(Nl,100,TR);

%%% DRAWING OF RAW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WhichTS2Draw=1;

YYY=ts(WhichTS2Draw,:);

% YYY = conv(YYY,hrf);
% YYY = YYY(1:Nl);

[Xx,Xy] = DrawMeSpectrum(X,2);
[Yx,Yy] = DrawMeSpectrum(YYY,2);

rawfh=figure('position',[50,500,600,900]);
subplot(5,1,1)
hold on; box on; grid on;
title('Resting-State Simulated Signal')
plot(YYY)
ylabel('a.u.')
xlabel(['Time (sampled on ' num2str(TR) 's)'])
subplot(5,1,2)
hold on; box on; grid on; 
title(['Design -- Design Efficiency: ' num2str(design_eff)])
plot(Xraw); plot(DD);
ylabel(['Time (sampled on ' num2str(TR) 's)'])
xlabel('Time')
subplot(5,1,3)
hold on; box on; grid on; 
title('Autocorrelation of Signal')
YYYac=autocorr(YYY,100); 
plot(YYYac)
ylabel('Autocorrelation Coefficients')
xlabel('Lags')
xlim([0 100])
subplot(5,1,4)
hold on; box on; grid on; 
title('Spectrum of the design')
plot(Xx,Xy)
ylabel('Magnitude')
xlabel('Frequency')
subplot(5,1,5)
hold on; box on; grid on; 
title('Spectrum of the Signal')
plot(Yx,Yy)
ylabel('Magnitude')
xlabel('Frequency')


%%% DRAWING OF Deterended %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

YYYd = HPFfilter*YYY';
Xd = HPFfilter*X;

[Xdx,Xdy] = DrawMeSpectrum(Xd,2);
[Ydx,Ydy] = DrawMeSpectrum(YYYd,2);

dfh=figure('position',[50,500,600,900]);
subplot(5,1,1)
hold on; box on; grid on;
title('Resting-State Simulated Signal -- DETERENDED')
plot(YYYd)
ylabel('a.u.')
xlabel(['Time (sampled on ' num2str(TR) 's)'])
subplot(5,1,2)
hold on; box on; grid on; 
title(['Design -- Design Efficiency: ' num2str(design_eff)])
plot(Xd); %plot(sf);
ylabel(['Time (sampled on ' num2str(TR) 's)'])
xlabel('Time')
subplot(5,1,3)
hold on; box on; grid on; 
title('Autocorrelation of Signal -- DETERENDED ')
YYYdac=autocorr(YYYd,100); 
plot(YYYdac)
ylabel('Autocorrelation Coefficients')
xlabel('Lags')
xlim([0 100])
subplot(5,1,4)
hold on; box on; grid on; 
title('Spectrum of the design -- DETERENDED')
plot(Xdx,Xdy)
ylabel('Magnitude')
xlabel('Frequency')
subplot(5,1,5)
hold on; box on; grid on; 
title('Spectrum of the Signal -- DETERENDED')
plot(Ydx,Ydy)
ylabel('Magnitude')
xlabel('Frequency')

%%% DRAWING OF PREW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
YYYd = YYYd - mean(YYYd);
Xd   = Xd - mean(Xd); 

% Get the residuals
Bhatd=pinv(Xd)*YYYd; % [b,c,d] = glmfit(Xdpw,YYYdpw);
YdHat=(Xd*Bhatd);
dRES=YYYd-YdHat; % sanitycheck: d.resid(1:5)

% Get the ACF of residual, regularised with a single Tukey of M
MTukey      = round(2*sqrt(Nl));
acfdRES     = AC_fft(dRES,Nl);
acfdRES_tt  = [1 TukeyTaperMe(acfdRES(2:end),Nl-1,MTukey)];
acfdRES_COV = toeplitz(acfdRES_tt); 

% ttfh=figure; 
% hold on; grid on; box on; 
% title(['Single Tukey Tapering of the ACF of residuals on M=' num2str(MTukey)])
% plot(acfdRES_tt)
% plot(acfdRES)
% ylabel('ACF'); xlabel('Lags')
% xlim([0 round(Nl/4)])

% Get the covariance matrix of residuals using Y-W AR1
%  AR1rho = AR_YW(dRES,1);
%  acfdRES_COV = full(spm_Q(AR1rho,Nl));

% prewhiten the signal & design with the covariance of residuals
YYYdpw = pinv(sqrtm(acfdRES_COV))*YYYd;
Xdpw = pinv(sqrtm(acfdRES_COV))*Xd;

% Recalculate the beta and the residuals of the prewhitened model 
Bhatdpw = pinv(Xdpw)*YYYdpw; % \Beta = X^{+}Y -- sanitycheck: [b,c,d] = glmfit(Xdpw,YYYdpw);
YdpwHat = Xdpw*Bhatdpw; %
dpwRES  = YYYdpw-YdpwHat; % sanitycheck: d.resid(1:5)

% Get the spectrums of signal
[YYYdx,YYYdy] = DrawMeSpectrum(YYYd,2);
[YYYdpwx,YYYdpwy] = DrawMeSpectrum(YYYdpw,2);
% Get the spectrums of residuals
[RESddx,RESddy] = DrawMeSpectrum(dRES,2);
[RESdpwdx,RESdpwdy] = DrawMeSpectrum(dpwRES,2);
% Get the spectrums of design
[Xdpwx,Xdpwy] = DrawMeSpectrum(Xdpw,2);

%DRAW

dfhpw=figure('position',[50,500,600,900]);
subplot(5,1,1)
hold on; box on; grid on;
title('Residual of the fit')
plot(dRES)
ylabel('a.u.')
xlabel(['Time (sampled on ' num2str(TR) 's)'])

subplot(5,1,2)
hold on; box on; grid on; 
title('Autocorrelation of Residual')
dRES_acf   = autocorr(dRES,100); 
dpwRES_acf = autocorr(dpwRES,100); 
plot(dRES_acf)
plot(dpwRES_acf)
ylabel('Autocorrelation Coefficients')
xlabel('Lags')
xlim([0 100])
legend({'Naive','Prewhitened'})

subplot(5,1,3)
hold on; box on; grid on; 
title(['Spectrum of Residuals'])
plot(RESddx,RESddy)
plot(RESdpwdx,RESdpwdy)
ylabel('Magnitude')
xlabel('Frequency')
legend({'Naive','Prewhitened'})

subplot(5,1,4)
hold on; box on; grid on; 
title('Spectrum of Design -- PREWHITENED')
plot(Xdx,Xdy)
plot(Xdpwx,Xdpwy)
ylabel('Magnitude')
xlabel('Frequency')
legend({'Naive','Prewhitened'})

% subplot(5,1,5)
% hold on; box on; grid on; 
% title('Spectrum of the Signal -- Deterend & PWed (ACF) ')
% plot(Ydpwdx,Ydpwdy)
% ylabel('Magnitude')
% xlabel('Frequency')


%%% SAVE FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(rawfh,'color','w'); export_fig(rawfh,'RFig/YXB.png')
set(dfh,'color','w'); export_fig(dfh,'RFig/YdXdB.png')
set(dfhpw,'color','w'); export_fig(dfhpw,'RFig/YpwXpwB.png')
set(ttfh,'color','w'); export_fig(ttfh,'RFig/ttacf.png')

% 
% frametimes=(0:Nt-1)*TR;
% slicetimes=[0.14 0.98 0.26 1.10 0.38 1.22 0.50 1.34 0.62 1.46 0.74 1.58 0.86];
% 
% eventid=kron(ones(10,1),[1; 2]);
% eventimes=(0:19)'*18+9;
% duration=ones(20,1)*9;
% height=ones(20,1);
% events=[eventid eventimes duration height];
% 
% X_cache=fmridesign(frametimes,slicetimes,events,[],hrf_parameters);