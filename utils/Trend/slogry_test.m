clear; 
%warning('off','all')

% COHORT = 'HCP'; 
% COHORTDIR = '/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/HCP/';
% pwdmethod = 'ACFadjT1S0'; %ACF AR-YW AR-W ARMAHR
% Mord      = 15; 
% lFWHM     = 0;
% SubID     = '101915';
% SesID     = 'REST1_LR'; 
% TR        = 0.720;
% EDtype    = 'ER'; % boxcar
% %TempTreMethod = 'spline'; 
% %NumTmpTrend   = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COHORT = 'ROCKLAND'; 
COHORTDIR = '/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/ROCKLAND/';
pwdmethod = 'ACF'; %ACF AR-YW AR-W ARMAHR
Mord      = 30; 
lFWHM     = 0;
SubID     = 'A00028150';
SesID     = 'DS2'; 
TR        = 0.645;
EDtype    = 'ER'; % boxcar
taskSNR   = 0.2;
icaclean  = 0;
%TempTreMethod = 'spline'; 
%NumTmpTrend   = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ACF:                  0.9917    0.9879    0.9879
%                       0.9889    0.9760    0.9762
%
% gACFadjxACFadj2tT1S0: 0.9708    0.9890    0.9887    (5mm)
%                       0.9028    0.9947    0.9934

TempTreMethod = 'hpf'; 
NumTmpTrend   = [];

%/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/R_ROCKLAND/R.PW/T900/ROCKLAND_645_900_ACFadjT0S0_AR--2_MA-0_FWHM5_hpfk_gsr0_aroma0/A00028150_DS

% What is flowing in from the cluster:
disp('From the cluster ======================')
disp(['SubID: ' SubID])
disp(['SesID: ' SesID])
disp(['TR: ' num2str(TR)])
disp(['ARmethod: ' pwdmethod])
disp(['AR order:' num2str(Mord)])
disp(['lFWHM: ' num2str(lFWHM)])
%disp(['COHORT directory:' COHORTDIR])




tic;
SaveImagesFlag      = 1; 
SaveMatFileFlag     = 1; 
DoDetrendingPrior   = 0; 
MParamNum           = 24;
gsrflag             = 0;
TempDerv            = 0;

disp('=======================================')

PATH2AUX='/Users/sorooshafyouni/Home/GitClone/FILM2';
addpath([PATH2AUX '/utils/Trend'])
addpath('/Users/sorooshafyouni/Home/matlab/spm12')
addpath([PATH2AUX '/utils/AR_YW'])
addpath([PATH2AUX '/utils/ARMA_HR'])
addpath([PATH2AUX '/mis'])


disp('=====SET UP PATHS =============================')
%Raw Images (MMP feat output)
Path2ImgRaw = [COHORTDIR '/R_mpp/sub-' SubID '/ses-' SesID];
if strcmpi(COHORT,'ROCKLAND')
    Path2ImgDir = [Path2ImgRaw '/sub-' SubID '_ses-' SesID '_task-rest_acq-' num2str(TR*1000) '_bold_mpp'];
elseif any(strcmpi(COHORT,{'Beijing','Cambridge'}))
    Path2ImgDir = [Path2ImgRaw '/rest_mpp'];
end

%Path2ImgRaw=[PATH2AUX '/ExampleData/R.mpp'];
Path2ImgDir = ['/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/ROCKLAND/sub-' SubID '/ses-' SesID '/sub-' SubID '_ses-' SesID '_task-rest_acq-645_bold_mpp'];

%Path2ImgDir = ['/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/' COHORT '/' SubID '_3T_rfMRI_' SesID '_mpp/'];

fwhmlab='';
% if lFWHM
%     fwhmlab=['_fwhm' num2str(lFWHM)];
% end

if ~icaclean || icaclean == -1
    icalab = 'off';
    Path2Img    = [Path2ImgDir '/prefiltered_func_data_bet' fwhmlab '.nii.gz'];
elseif icaclean==1
    icalab = 'nonaggr';
    Path2Img    = [Path2ImgDir '/ica-aroma_hpf' fwhmlab '/denoised_func_data_nonaggr.nii.gz'];
elseif icaclean==2
    icalab = 'aggr';
    Path2Img    = [Path2ImgDir '/ica-aroma_hpf' fwhmlab '/denoised_func_data_aggr.nii.gz'];
end

path2mask = [Path2ImgDir '/mask.nii.gz']; 
Path2MC   = [Path2ImgDir '/prefiltered_func_data_mcf.par'];

SEGmaskinEPI = [Path2ImgDir '/seg/sub-' SubID '_ses-' SesID '_T1w_func_seg.nii.gz'];
WMseg        = [Path2ImgDir '/seg/sub-' SubID '_ses-' SesID '_T1w_func_seg_2.nii.gz'];

disp(['Image: ' Path2Img])
disp(['Motion params: ' Path2MC])

% Directory 2 save the results
Path2ImgResults=[PATH2AUX '/ExampleData/R.mpp/RNullfMRI_' SubID '_' SesID];
if ~exist(Path2ImgResults, 'dir')
	mkdir(Path2ImgResults)
	disp(['The directory: ' Path2ImgResults ' did not exists. I made one. '])
end

disp(['Output stuff: ' Path2ImgResults])

%%% Read The Data %%%%%%%%%%%%%%%%%%%%%%%%
disp('=====LOAD THE IMAGE ===========================')

[Y,InputImgStat] = CleanNIFTI_spm(Path2Img,'demean','norm',100);
Y = Y';
T = size(Y,1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y_ACL                  = ACLImage(Y'); 

trds                   = sgolayfilt(Y,3,241);
% [~,PCs,~,~,trds_var]   = pca(trds-mean(trds));
% trds_var               = cumsum(trds_var); 
% PCs                    = PCs(:,(trds_var<95));
% [~,~,YY]               = myOLS(Y,PCs);

YY                     = Y - trds; 

YY_ACL                 = ACLImage(YY'); 

%%%%%%%%%%%%%% HPF %%%%%%%%%%%%%%%%%%%%%%%%%%%

hp_ff  = hp_fsl(T,100,TR);    
YY0     = hp_ff*Y;  % high pass filter the data    
YY0_ACL = ACLImage(YY0'); 














