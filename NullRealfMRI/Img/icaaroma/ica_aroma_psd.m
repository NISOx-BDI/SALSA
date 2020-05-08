clear; 

SubID     = 'A00008399';
SesID     = 'DS2'; 
TR        = 0.645;
MParamNum = 24; 
disp('=======================================')

PATH2AUX='/Users/sorooshafyouni/Home/GitClone/FILM2';
addpath([PATH2AUX '/utils/Trend'])
addpath('/Users/sorooshafyouni/Home/matlab/spm12')
addpath([PATH2AUX '/utils/AR_YW'])
addpath([PATH2AUX '/utils/ARMA_HR'])
addpath([PATH2AUX '/utils/Spectral'])
addpath([PATH2AUX '/mis'])


disp('=====SET UP PATHS =============================')
%Raw Images (MMP feat output)
Path2ImgRaw=[PATH2AUX '/ExampleData/R.mpp'];
%Path2ImgDir = ['/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/R_test/sub-' SubID '_ses-' SesID '_task-rest_acq-645_bold_mpp'];
Path2ImgDir = ['/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/ROCKLAND/sub-' SubID '/ses-' SesID '/sub-' SubID '_ses-' SesID '_task-rest_acq-645_bold_mpp'];
Path2Img_ica_aroma    = [Path2ImgDir '/ica-aroma'];
Path2Img_melodic    = [Path2Img_ica_aroma '/melodic.ica'];
Path2MC  = [Path2ImgDir '/prefiltered_func_data_mcf.par'];

X   = []; 
MCp = [];
if MParamNum     == 6
    MCp = load(Path2MC);
elseif MParamNum == 12
    MCp = load(Path2MC);
    MCp = [MCp,MCp.^2]; % 12 parameter motion 
elseif MParamNum == 24
    o6MCp   = load(Path2MC);  % 6 orig param
    so6MCp  = o6MCp.^2;       % 6 square orig
    do6MCp  = diff(o6MCp);    % 6 diff
    do6MCp  = [ do6MCp; zeros(1,size(do6MCp,2)) ];
    sd6MCp  = do6MCp.^2;      % 6 square of diff
    %sd6MCp  = [ sd6MCp; zeros(1,size(sd6MCp,2)) ];
    MCp     = [o6MCp,so6MCp,do6MCp,sd6MCp];
end
disp(['Number of motion parameter: ' num2str(MParamNum)])

X = [X,MCp];

% noise ICs
nic = load([Path2Img_ica_aroma '/classified_motion_ICs.txt']);

% time courses
ts = load([Path2Img_melodic '/melodic_mix']);
[T,numIC] = size(ts); 
% Regress out the motions
pinvX           = pinv(X); 
ResidFormingMat = eye(T)-X*pinvX; % residual forming matrix 
ts          = ResidFormingMat*ts;
% high pass filter
hp_ff = hp_fsl(T,100,TR);    
ts    = hp_ff*ts;    % high pass filter the design


ts = ts-mean(ts); 
ts = ts./std(ts); 

% noise time series
nts = ts(:,nic);

% signal time series
sts = ts; 
sts(:,nic) = []; 

%% draw the spd
% [sts_freq,sts_spd] = DrawMeSpectrum(sts,TR); 
% [nts_freq,nts_spd] = DrawMeSpectrum(nts,TR); 
% 
% figure; hold on; 
% plot(sts_freq,mean(sts_spd,2),'r')
% plot(nts_freq,mean(nts_spd,2),'b')
% 
% %% prewhiten
% ntspwd = PrewhitenSPD(nts,30,1/TR); 
% stspwd = PrewhitenSPD(sts,30,1/TR); 
% 
% [sts_freq,sts_spd] = DrawMeSpectrum(stspwd,TR); 
% [nts_freq,nts_spd] = DrawMeSpectrum(ntspwd,TR); 
% 
% figure; hold on; 
% plot(sts_freq,mean(sts_spd,2),'r')
% plot(nts_freq,mean(nts_spd,2),'b')
