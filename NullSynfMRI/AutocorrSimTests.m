clear

SubID     = 'A00027167';
SesID     = 'DS2'; 
TR        = 0.645;
disp('=======================================')
PATH2AUX='/Users/sorooshafyouni/Home/GitClone/FILM2';
addpath([PATH2AUX '/utils/Trend'])
addpath([PATH2AUX '/mis'])

addpath('/Users/sorooshafyouni/Home/matlab/spm12')
%12487 -- 7061
disp('=====SET UP PATHS =============================')
%Raw Images (MMP feat output)
Path2ImgRaw=[PATH2AUX '/ExampleData/R.mpp'];
%Path2ImgDir = [Path2ImgRaw '/sub-' SubID '/ses-' SesID '/sub-' SubID '_ses-' SesID '_task-rest_acq-' num2str(TR*1000) '_bold.mpp'];
Path2ImgDir = ['/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/R_test/sub-' SubID '_ses-' SesID '_task-rest_acq-645_bold_mpp'];
Path2Img    = [Path2ImgDir '/prefiltered_func_data_bet.nii'];


[Y,S] = CleanNIFTI_spm(Path2Img,'demean');
T     = size(Y,2);
Vorig = size(Y,1); 

if size(Y,1)~=T; Y = Y'; end; %TxV

pp = 1+fix(TR.*T/150);
dY = multpolyfit(repmat(1:T,Vorig,1)',Y,T,pp)';
dY = dY - mean(dY); 

acf = AC_fft(dY,T);
acl = sum(acf.^2,2);

[Xp,Yp,specflat] = DrawMeSpectrum(dY,TR);

xx = dY(:,randi(Vorig));
ac = autocorr(xx,T-1);

[InvMat,spd] = QuickAutoCorr(ac(1:50),T);
spd
sxx = InvMat*randn(T,1); 
sac = autocorr(sxx,T-1);

figure; hold on; 
plot(sac(1:200))
plot(ac(1:200))

%scatter(1./acl,specflat)
