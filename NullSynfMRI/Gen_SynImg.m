clear; warning('off','all')

pwdmethod = 'ACF'; %ACF AR-YW AR-W ARMAHR
Mord      = 30; 
lFWHM     = 0;
SubID     = 'A00028858';
SesID     = 'DS2'; 
TR        = 0.645;
T         = 500;
SimMord   = 50;

% What is flowing in from the cluster:
disp('From the cluster ======================')
disp(['SubID: ' SubID])
disp(['SesID: ' SesID])
disp(['TR: ' num2str(TR)])
disp(['ARmethod: ' pwdmethod])
disp(['AR order:' num2str(Mord)])
disp(['lFWHM: ' num2str(lFWHM)])
%disp(['COHORT directory:' COHORTDIR])
disp('=======================================')

PATH2AUX='/Users/sorooshafyouni/Home/GitClone/FILM2';
addpath([PATH2AUX '/utils/Trend'])
addpath('/Users/sorooshafyouni/Home/matlab/spm12')
addpath([PATH2AUX '/utils/AR_YW'])
addpath([PATH2AUX '/utils/ARMA_HR'])
addpath([PATH2AUX '/mis'])

disp('=====SET UP PATHS =============================')
%Raw Images (MMP feat output)
Path2ImgRaw=[PATH2AUX '/ExampleData/R.mpp'];
%Path2ImgDir = [Path2ImgRaw '/sub-' SubID '/ses-' SesID '/sub-' SubID '_ses-' SesID '_task-rest_acq-' num2str(TR*1000) '_bold.mpp'];
Path2ImgDir = ['/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/ROCKLAND/sub-' SubID '/ses-' SesID '/sub-' SubID '_ses-' SesID '_task-rest_acq-645_bold_mpp'];

if ~lFWHM
    Path2Img = [Path2ImgDir '/prefiltered_func_data_gm_bet.nii.gz'];
else
    Path2Img = [Path2ImgDir '/prefiltered_func_data_bet_FWHM' num2str(lFWHM) '.nii'];
end

Path2MC  = [Path2ImgDir '/prefiltered_func_data_mcf.par'];

disp(['Image: ' Path2Img])
disp(['Motion params: ' Path2MC])

% Directory 2 save the results
Path2ImgResults=[PATH2AUX '/ExampleData/R.mpp/RNullfMRI_' SubID '_' SesID];

if ~exist(Path2ImgResults, 'dir')
	mkdir(Path2ImgResults)
	disp(['The directory: ' Path2ImgResults 'did not exists. I made one. '])
end

disp(['Output stuff: ' Path2ImgResults])

SaveMatFileFlag     = 1;
SaveImagesFlag      = 1; 
MParamNum           = 24; 

%%% Read The Data %%%%%%%%%%%%%%%%%%%%%%%%
disp('=====LOAD THE IMAGE ===========================')
[Y,InputImgStat]=CleanNIFTI_spm(Path2Img,'demean');
Ti = InputImgStat.CleanedDim(2);
TR = InputImgStat.voxelsize(4);
Vorig = InputImgStat.CleanedDim(1);
V = Vorig;

if size(Y,1)~=Ti; Y = Y'; end; %TixV
%%% DETREND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NumTmpTrend = 1+fix(TR.*Ti/150);
dY = multpolyfit(repmat(1:Ti,Vorig,1),Y',Ti,NumTmpTrend)';
dY = dY - mean(dY); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DESIGN MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('++++++++++++ Construct a design matrix')
%%% Generate a Design Matrix --------------------------------
EDtype = 'boxcar'; 
BCl = 20;
EDX = GenerateED(BCl,T,TR); 
EDX = EDX - mean(EDX); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Motion Params %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Motion parameters 
o6MCp   = load(Path2MC);  % 6 orig param
so6MCp  = o6MCp.^2;       % 6 square orig
do6MCp  = diff(o6MCp);    % 6 diff
do6MCp  = [ do6MCp; zeros(1,size(do6MCp,2)) ];
sd6MCp  = do6MCp.^2;      % 6 square of diff
%sd6MCp  = [ sd6MCp; zeros(1,size(sd6MCp,2)) ];
MCp     = [o6MCp,so6MCp,do6MCp,sd6MCp];

disp(['Number of motion parameter: ' num2str(MParamNum)])

X           = MCp;
X           = X - mean(X); % demean everything 
%%% RESIDUALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pinvX            = pinv(X); 
ResidFormingMat0 = eye(Ti)-X*pinvX; % residual forming matrix 
dY               = ResidFormingMat0*dY;

pinvEDX           = pinv(EDX); 
ResidFormingMat  = eye(T)-EDX*pinvEDX; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BIAS REDUCTION OF AUTOREGRESSIVE MODELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ACFadjflag = 0; WrosleyFlag = 0; ACFflag = 0; ARMAHRflag = 0; MPparamNum = 0; 
if strcmpi(pwdmethod,'AR-W') %Worsely %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    WrosleyFlag = 1;     
    invM_biasred = ACFBiasAdj(ResidFormingMat,T,Mord);    
elseif strcmpi(pwdmethod,'ACFadj') % Yule-Walker %%%%%%%%%%%%%%%%%%%%%%%%%%%
    ACFadjflag          = 1; 
    invM_biasred = ACFBiasAdj(ResidFormingMat,T,Mord); 
elseif strcmpi(pwdmethod,'ACF') % ACF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ACFflag         = 1;
elseif strcmpi(pwdmethod,'ARMAHR') % ARMAHR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ARMAHRflag      = 1; 
    MPparamNum      = 1;  % the MA order 
    ARParamARMA     = 50; % the higher fit in ARMA HR
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PREWHITEN THE RESIDULAS & ESTIMATE BIAS AND CPS %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% ACFs %%%%%%%%%%%%%%%
disp(['Calculate the autocorrelation coefficients.'])
[~,~,dYacov]  = AC_fft(dY,Ti); % Autocovariance 
dYacov        = dYacov'; %TxV
dYacorr       = dYacov./sum(abs(dY).^2); % Autocorrelation
ACL           = sum(dYacorr.^2); % Autocorrelation Length

%% RAM SAVE %%
clear dY Y
%%%%%%%%%%%%%%

%%% Preallocate memory
Bhat_PW    = zeros(V,1);
SE_PW      = zeros(V,1); 
tVALUE_PW  = zeros(V,1);
CPSstat_PW = zeros(V,1); 
CPZ_PW     = zeros(V,1);
dpwRES     = zeros(T,V); 
sY         = zeros(T,V);

for vi = 1:V
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Simulate the time series %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
%     %Generate a time series -------------------------   
    [ACMatDecomp,cholord]  = QuickAutoCorr(dYacorr(1:SimMord,vi),T); 
    
    if cholord<Mord; disp('cholesky failed on minimum AR order requested.'); end; 
    
    z            = randn(T,1);
    z            = z-mean(z); 
    sY_tmp       = ACMatDecomp*z;
    
    sY(:,vi)     = sY_tmp;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Prewhiten the residuals %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

    residY      	    = ResidFormingMat*sY_tmp;
    residY     	 	    = residY-mean(residY); 
    [dRESacorr,~,dRESacov]  = AC_fft(residY',T);
    dRESacorr = dRESacorr';
    dRESacov  = dRESacov';
    
    %dRESacorr  = autocorr(residY,T-1); % pretty fast on single time series, autocorr doesn't exist in Octave -- WTF!
    %dRESacov   = dRESacorr.*sum(abs(residY).^2);
    
    if ACFadjflag % Yule-Walker %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~mod(vi,5000); disp([ pwdmethod ' ::: on voxel ' num2str(vi)]); end;        
%         AR_YW -------------------------------------------
%         ac_tmp          = dRESacorr;    
%         R_tmp           = toeplitz(ac_tmp(1:Mord));
%         r_tmp           = ac_tmp(2:Mord+1);
%         %YWARparam_tmp   = pinv(R_tmp)*r_tmp;
%         YWARparam_tmp   = R_tmp\r_tmp;
%         % ------------------------------------------------
%         ACMat           = full(spm_Q(YWARparam_tmp,T));
%         [sqrtmVhalf,spdflag] = CholWhiten(ACMat);
%         %invACMat        = inv(ACMat); % pinv is damn slow!
%         %sqrtmVhalf      = chol(invACMat);

        [sqrtmVhalf,spdflag] = ACF_ResPWm(dRESacov,Mord,invM_biasred,1);


    elseif ACFflag % ACF - Tukey tapered %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~mod(vi,5000); disp([ pwdmethod ' ::: on voxel ' num2str(vi)]); end;        
%         acfdRES     = dRESacorr; 
%         %Mord        = 2*round(sqrt(T));
%         acfdRES_tt  = [1 TukeyTaperMe(acfdRES(2:end),T-1,Mord)];
%         ACMat       = toeplitz(acfdRES_tt);  
%         [sqrtmVhalf,spdflag] = CholWhiten(ACMat);
        
        [sqrtmVhalf,spdflag] = ACF_ResPWm(dRESacov,Mord,[],1);

    elseif WrosleyFlag % Worsely %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~mod(vi,5000); disp([ pwdmethod ' ::: on voxel ' num2str(vi)]); end;        
%         dRESac_adj      = (invM_biasred*dRESacov(1:Mord+1));
%         dRESac_adj      = dRESac_adj./dRESac_adj(1); % make a auto*correlation*
%         [Ainvt,posdef]          = chol(toeplitz(dRESac_adj)); 
%         p1                      = size(Ainvt,1); 
%         A                       = inv(Ainvt'); 
%         sqrtmVhalf              = toeplitz([A(p1,p1:-1:1) zeros(1,T-p1)],zeros(1,T)); 
%         sqrtmVhalf(1:p1,1:p1)   = A;

           [sqrtmVhalf,spdflag] = AR_ResPWm(dRESacov,Mord,invM_biasred);

    elseif ARMAHRflag % ARMA HR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~mod(vi,5000); disp([ pwdmethod ' ::: on voxel ' num2str(vi)]); end;        
  
        %AR_YW -------------------------------------------
        ac_tmp              = dRESacorr;    
        R_tmp               = toeplitz(ac_tmp(1:ARParamARMA));
        r_tmp               = ac_tmp(2:ARParamARMA+1);
        %YWARparam_tmp       = pinv(R_tmp)*r_tmp %
        YWARparam_tmp       = R_tmp\r_tmp;
        % ------------------------------------------------
        
        [arParam,maParam]    = ARMA_HR_ACF(residY,YWARparam_tmp',T,Mord,MPparamNum);
        ACMat                = ARMACovMat([arParam,maParam],T,Mord,MPparamNum);
        [sqrtmVhalf,spdflag] = CholWhiten(ACMat);
    end
    
    % Make the X & Y whitened %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ystar_YW = sqrtmVhalf*sY_tmp;
    Xstar_YW = sqrtmVhalf*EDX;
    
    % Fit a model to the prewhitened system  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    Xstar_YW                                      = [ones(T,1), Xstar_YW]; % add intercept
    [Bhat_PW_S_tmp,~,dpwRES_tmp,Stat_PW_SE_T_tmp] = myOLS(Ystar_YW,Xstar_YW,[0 1]);
    Bhat_PW(vi)    = Bhat_PW_S_tmp;         
    SE_PW(vi)      = Stat_PW_SE_T_tmp.se;
    tVALUE_PW(vi)  = Stat_PW_SE_T_tmp.tval;

    dpwRES(:,vi)       = dpwRES_tmp;
    
    [~,CPSstat_PW(vi),CPZ_PW(vi)] = CPSUnivar(dpwRES_tmp,Xstar_YW);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Naive Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('++++++++++++Fit data to the Naive model.')
X0                                          = [ones(T,1),EDX];
[Bhat_Naive,~,resNaive,Stat_Naive_SE_tmp]   = myOLS(sY,X0,[0 1]);
SE_Naive                                    = Stat_Naive_SE_tmp.se;
tVALUE_Naive                                = Stat_Naive_SE_tmp.tval;
[~,CPSstat_Naive,CPZ_Naive]                 = CPSUnivar(resNaive,X0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECTRUM OF THE RESIDUALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('++++++++++++Calculate the spectrum of the residuals.')
[dpwRESXp,dpwRESYp,fs_PW] = DrawMeSpectrum(dpwRES,TR,0);
dpwRESYp            = mean(dpwRESYp,2); % average across voxels

clear dpwRES

[resNaiveSXp,resNaiveYp,fs_Naive] = DrawMeSpectrum(resNaive,TR,0);
resNaiveYp               = mean(resNaiveYp,2); % average across voxels


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE THE RESULTS AS AN IMAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SaveImagesFlag
    % 3D IMAGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    VariableList = {'SE_PW','SE_Naive',...
		'Bhat_Naive','Bhat_PW',...
		'tVALUE_Naive','tVALUE_PW',...
 		'fs_PW','fs_Naive'};
    OutputImgStat            = InputImgStat.spmV(1);
    OutputImgStat.Removables = InputImgStat.Removables;

    for vname = VariableList
        tmpvar                   = eval(vname{1});
        OutputImgStat.fname      = [Path2ImgResults '/Sim'  num2str(SIDX)  '_ED' EDtype '_' num2str(BCl) '_' pwdmethod '_AR' num2str(Mord) '_MA' num2str(MPparamNum)  '_' vname{1} '.nii'];

        CleanNIFTI_spm(tmpvar,'ImgInfo',OutputImgStat);
	system(['gzip ' OutputImgStat.fname])
    end
end

% MAT FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SaveMatFileFlag
    GLM.df = Stat_Naive_SE_tmp.df;
    GLM.X  = X0;
    GLM.C  = [0 1];
    GLM.EDtype = EDtype;
    GLM.EDFreq = BCl;

    SPEC.X_RES   = resNaiveSXp;
    SPEC.Y_RES   = resNaiveYp;
    SPEC.X_pwRES = dpwRESXp;
    SPEC.Y_pwRES = dpwRESYp;

    PW.dt     = 'poly';
    PW.dtl    = NumTmpTrend;
    PW.pwmeth = pwdmethod;
    PW.fwhm   = lFWHM;
    PW.MAp    = MPparamNum;
    PW.ARp    = Mord;

    MatFileName = [Path2ImgResults '/Sim'  num2str(SIDX)  '_ED' EDtype '_' num2str(BCl) '_' pwdmethod '_AR' num2str(Mord) '_MA' num2str(MPparamNum) '_FWHM' num2str(lFWHM) '.mat'];
    save(MatFileName,'GLM','SPEC','PW')
end

disp('--DONE--')


