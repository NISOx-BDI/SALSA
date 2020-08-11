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
pwdmethod = 'ACFadjT1S0'; %ACF AR-YW AR-W ARMAHR
Mord      = 10; 
lFWHM     = 0;
SubID     = 'A00028150';
SesID     = 'DS2'; 
TR        = 0.645;
EDtype    = 'ER'; % boxcar
%TempTreMethod = 'spline'; 
%NumTmpTrend   = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


TempTreMethod = 'hpf'; 
NumTmpTrend   = [];

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
if lFWHM
    fwhmlab=['_fwhm' num2str(lFWHM)];
end

Path2Img    = [Path2ImgDir '/prefiltered_func_data_bet' fwhmlab '.nii.gz'];

path2mask = [Path2ImgDir '/mask.nii.gz']; 
Path2MC   = [Path2ImgDir '/prefiltered_func_data_mcf.par'];

%%% Read The Data %%%%%%%%%%%%%%%%%%%%%%%%
disp('=====LOAD THE IMAGE ===========================')

[Y,InputImgStat]=CleanNIFTI_spm(Path2Img,'demean','norm',100);

%---------------------------------------------

T     = InputImgStat.CleanedDim(2);
TR    = InputImgStat.voxelsize(4);
Vorig = InputImgStat.CleanedDim(1);
V     = Vorig;

if size(Y,1)~=T; Y = Y'; end; %TxV

TempTrend = []; K = []; 
if ~exist('NumTmpTrend','var'); NumTmpTrend=[]; end;
if any(strcmpi(TempTreMethod,{'dct','spline','poly'}))
    [TempTrend,NumTmpTrend]   = GenerateTemporalTrends(T,TR,TempTreMethod,NumTmpTrend); % DC + Poly trends + Spline trends 
    TempTrend   = TempTrend(:,2:end); % we add a column of one later.
elseif strcmpi(TempTreMethod,{'hpf'})
    if isempty(NumTmpTrend) || ~exist('NumTmpTrend','var'); NumTmpTrend=100; end; 
    hp_ff = hp_fsl(T,NumTmpTrend,TR);    
    Y     = hp_ff*Y;  % high pass filter the data
end
disp(['Detrending: ' TempTreMethod ',param: ' num2str(NumTmpTrend)])


% Save the image
% CleanNIFTI_spm(tmpvar,'ImgInfo',InputImgStat.spmV,'DestDir',fname,'removables',InputImgStat.Removables);
fname = [Path2ImgDir '/filtered_' TempTreMethod '_func_data_bet' fwhmlab '.nii']; 
[Y,InputImgStat]=CleanNIFTI_spm(Y','ImgInfo',InputImgStat.spmV,'DestDir',fname,'removables',InputImgStat.Removables);
