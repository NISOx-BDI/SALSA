%clear

warning('off','all')


% ---------------- TEST ----------------------------
%pwdmethod  = 'ACF'; %ACF AR-YW AR-W ARMAHR
%Mord       = 30; 
%MPparamNum = 0;
%TempTreMethod = 'dct'
%SubID     = 'A00029304';
%SesID     = 'DS2';
%lFWHM     = 0;
%TR        = 0.645;
%COHORTDIR = '/well/nichols/users/scf915/ROCKLAND';
%Path2ImgResults = [COHORTDIR '/R.PW/TEST' pwdmethod '_AR-' num2str(Mord) '_MA-' num2str(MPparamNum)  ]




% What is flowing in from the cluster:
disp('From the cluster ======================')
disp(['SubID: ' SubID])
disp(['SesID: ' SesID])
disp(['TR: ' num2str(TR)])
disp(['ARmethod: ' pwdmethod])
disp(['AR order:' num2str(Mord)])
disp(['MA order: ' num2str(MPparamNum)])
disp(['lFWHM: ' num2str(lFWHM)])
disp(['Detrending: ' TempTreMethod])
disp(['COHORT directory:' COHORTDIR])
disp(['Parth 2 Results: ' Path2ImgResults ])
disp('=======================================')

PATH2AUX='~/bin/FILM2';
addpath([PATH2AUX '/utils/Trend'])
addpath([PATH2AUX '/utils/AR_YW'])
addpath([PATH2AUX '/utils/ARMA_HR'])
addpath([PATH2AUX '/mis'])
addpath (fullfile ('/users/nichols/scf915', 'spm12-r7771'));
%addpath('/well/nichols/users/scf915/externals/spm12')

disp('=====SET UP PATHS =============================')
%Raw Images (MMP feat output)
Path2ImgRaw = [COHORTDIR '/R_mpp'];
Path2ImgDir = [Path2ImgRaw '/sub-' SubID '/ses-' SesID '/sub-' SubID '_ses-' SesID '_task-rest_acq-' num2str(TR*1000) '_bold_mpp'];

if ~lFWHM
    Path2Img    = [Path2ImgDir '/prefiltered_func_data_bet.nii'];
else
    Path2Img    = [Path2ImgDir '/prefiltered_func_data_bet_FWHM' num2str(lFWHM) '.nii'];
end
  
Path2MC     = [Path2ImgDir '/prefiltered_func_data_mcf.par'];

disp(['Image: ' Path2Img])
disp(['Motion params: ' Path2MC])

% Directory 2 save the results
Path2ImgResults=[Path2ImgResults '/' SubID '_' SesID];
if ~exist(Path2ImgResults, 'dir')
	mkdir(Path2ImgResults)
	disp(['The directory: ' Path2ImgResults ' did not exists. I made one. '])
end

disp(['Output stuff: ' Path2ImgResults])

SaveImagesFlag      = 1; 
SaveMatFileFlag     = 1; 
DoDetrendingPrior   = 0; 
MParamNum           = 24; 

%%% Read The Data %%%%%%%%%%%%%%%%%%%%%%%%
disp('=====LOAD THE IMAGE ===========================')
[Y,InputImgStat]=CleanNIFTI_spm(Path2Img,'demean');
T = InputImgStat.CleanedDim(2);
TR = InputImgStat.voxelsize(4);
Vorig = InputImgStat.CleanedDim(1);
V = Vorig;

if size(Y,1)~=T; Y = Y'; end; %TxV

%%% DETREND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DoDetrendingPrior
    pp = 1+fix(TR.*T/150);
    dY = multpolyfit(repmat(1:T,Vorig,1),Y,T,pp);
    dY = dY - mean(dY,2); 
else
    dY = Y; 
end

%%% SIMULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Order of AR which will be used to simulated stuff
% ARporder = 20;
% ARParam  = AR_YW(dYorig,ARporder);
% nRlz = 1000; 
% model = arima('Constant',0,'AR',ARParam,'Variance',1);
% dY = simulate(model,T,'NumPaths',nRlz)'; 
% dY = dY - mean(dY,2); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DESIGN MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('++++++++++++ Construct a design matrix')
%%% Generate a Design Matrix --------------------------------
EDtype = 'boxcar'; 
BCl = 20;
EDX = GenerateED(BCl,T,TR); 
EDX = EDX - mean(EDX); 

X   = EDX;
Xc  = 1; % where is the experimental design?

disp(['design updated, ' num2str(size(X,2))])
% Motion parameters ----------------------------------------
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
disp(['design updated, ' num2str(size(X,2))])

% Temporal trends ----------------------------------------
TempTrend = [];
if ~exist('NumTmpTrend','var'); NumTmpTrend=[]; end;
if any(strcmpi(TempTreMethod,{'dct','spline','poly'}))
    [TempTrend,NumTmpTrend]   = GenerateTemporalTrends(T,TR,TempTreMethod,NumTmpTrend); % DC + Poly trends + Spline trends 
    TempTrend   = TempTrend(:,2:end); % we add a column of one later.
elseif strcmpi(TempTreMethod,{'hpf'})
    if isempty(NumTmpTrend) || ~exist('NumTmpTrend','var'); NumTmpTrend=100; end; 
    hp_ff = hp_fsl(T,NumTmpTrend,TR);    
    X     = hp_ff*X;    % high pass filter the design
    dY    = hp_ff*dY;  % high pass filter the data
end
disp(['Detrending: ' TempTreMethod ',param: ' num2str(NumTmpTrend)])
%
X           = [X,TempTrend];
disp(['design updated, ' num2str(size(X,2))])

% Centre the design  ----------------------------------
X           = X - mean(X); % demean everything 

%%% RESIDUALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('++++++++++++Get the residuals using I-XX+')
pinvX           = pinv(X); 
ResidFormingMat = eye(T)-X*pinvX; % residual forming matrix 
residY          = ResidFormingMat*dY;

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
%%% FIT A MODEL TO THE ACd DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('++++++++++++Fit data to the Naive model.')
X0                                          = [ones(T,1),X];
glmcont                                     = zeros(1,size(X0,2));
glmcont(Xc+1)                               = 1;
[Bhat_Naive,~,resNaive,Stat_Naive_SE_tmp]   = myOLS(dY,X0);
SE_Naive                                    = Stat_Naive_SE_tmp.se;
tVALUE_Naive                                = Stat_Naive_SE_tmp.tval;
[~,CPSstat_Naive,CPZ_Naive]                 = CPSUnivar(resNaive,X0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AUTOCORR & AUTOCOV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% ACFs %%%%%%%%%%%%%%%
disp(['++++++++++++Calculate the autocorrelation coefficients.'])
%residY          = residY-mean(residY); 
residY          = residY-repmat(mean(residY),T,1);
[~,~,dRESacov]  = AC_fft(residY,T); % Autocovariance; VxT
dRESacov        = dRESacov'; %TxV
dRESacorr       = dRESacov./sum(abs(residY).^2); % Autocorrelation
ACL             = sum(dRESacorr.^2); % Autocorrelation Length

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PREWHITEN THE RESIDULAS & ESTIMATE BIAS AND CPS %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Preallocate memory
Bhat_PW    = zeros(V,1);
SE_PW      = zeros(V,1); 
tVALUE_PW  = zeros(V,1);
CPSstat_PW = zeros(V,1); 
CPZ_PW     = zeros(V,1);
dpwRES     = zeros(T,V); 
nonstationaryvox = [];
disp('++++++++++++Starts the voxel-wise prewhitening')

for vi = 1:V
    spdflag = 0;
    if ACFadjflag % ACF, Tapered, Adjusted %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~mod(vi,5000); disp([ pwdmethod ' ::: on voxel ' num2str(vi)]); end; 
        
        [sqrtmVhalf,spdflag] = ACF_ResPWm(dRESacov(:,vi),Mord,invM_biasred,1);
           
    elseif ACFflag % ACF - Tapered, Adjusted %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~mod(vi,5000); disp([ pwdmethod ' ::: on voxel ' num2str(vi)]); end; 
        
        [sqrtmVhalf,spdflag] = ACF_ResPWm(dRESacov(:,vi),Mord,[],1);
        
    elseif WrosleyFlag % Worsely %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~mod(vi,5000); disp([ pwdmethod ' ::: on voxel ' num2str(vi)]); end; 
        
        [sqrtmVhalf,spdflag] = AR_ResPWm(dRESacov(:,vi),Mord,invM_biasred);
        
    elseif ARMAHRflag % ARMA HR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~mod(vi,5000); disp(['ARMA-HR ::: on voxel ' num2str(vi)]); end; 
  
        %AR_YW -------------------------------------------
        ac_tmp              = dRESacorr(:,vi);    
        R_tmp               = toeplitz(ac_tmp(1:ARParamARMA));
        r_tmp               = ac_tmp(2:ARParamARMA+1);
        %YWARparam_tmp       = pinv(R_tmp)*r_tmp %
        YWARparam_tmp       = R_tmp\r_tmp;
        % ------------------------------------------------
        
        [arParam,maParam]    = ARMA_HR_ACF(residY(:,vi),YWARparam_tmp',T,Mord,MPparamNum);
        ACMat                = ARMACovMat([arParam,maParam],T,Mord,MPparamNum);
        [sqrtmVhalf,spdflag] = CholWhiten(ACMat);
    end

    if spdflag
        nonstationaryvox = [nonstationaryvox vi];
    end
    % Make the X & Y whitened %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ystar_YW = sqrtmVhalf*dY(:,vi);
    Xstar_YW = sqrtmVhalf*X;   
    %YvWY(vi) = corr(Ystar_YW,dY(:,vi)); % idon't know how useful that is.
    % Fit a model to the prewhitened system  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    Xstar_YW                                      = [ones(T,1), Xstar_YW]; % add intercept
    [Bhat_PW_S_tmp,~,dpwRES_tmp,Stat_PW_SE_T_tmp] = myOLS(Ystar_YW,Xstar_YW,glmcont);
    Bhat_PW(vi)    = Bhat_PW_S_tmp;         
    SE_PW(vi)      = Stat_PW_SE_T_tmp.se;
    tVALUE_PW(vi)  = Stat_PW_SE_T_tmp.tval;

    dpwRES(:,vi)       = dpwRES_tmp;
    
    [~,CPSstat_PW(vi),CPZ_PW(vi)] = CPSUnivar(dpwRES_tmp,Xstar_YW);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECTRUM OF THE RESIDUALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('++++++++++++Calculate the spectrum of the residuals.')
[dpwRESXp,dpwRESYp] = DrawMeSpectrum(dpwRES,TR,0);
dpwRESYp            = mean(dpwRESYp,2); % average across voxels

clear dpwRES

[resNaiveSXp,resNaiveYp] = DrawMeSpectrum(resNaive,TR,0);
resNaiveYp               = mean(resNaiveYp,2); % average across voxels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE THE RESULTS AS AN IMAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('++++++++++++Save the results.')

if SaveImagesFlag
    % 3D IMAGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    VariableList = {'Bhat_Naive','SE_Naive','tVALUE_Naive',...
        'Bhat_PW','SE_PW','tVALUE_PW',...
        'CPSstat_PW','CPZ_PW',...
        'CPZ_Naive','CPSstat_Naive',...
        'ACL'};
    OutputImgStat            = InputImgStat.spmV(1);
    OutputImgStat.Removables = InputImgStat.Removables;

    for vname = VariableList

        tmpvar                   = eval(vname{1});
        OutputImgStat.fname      = [Path2ImgResults '/sub-' SubID '_ses-' SesID '_ED' EDtype '_' num2str(BCl) '_' pwdmethod '_AR' num2str(Mord) '_MA' num2str(MPparamNum) '_FWHM' num2str(lFWHM) '_' TempTreMethod num2str(NumTmpTrend) '_' vname{1} '.nii'];

        CleanNIFTI_spm(tmpvar,'ImgInfo',OutputImgStat);
        system(['gzip ' OutputImgStat.fname]);
    end
end

% MAT FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SaveMatFileFlag
    GLM.df = Stat_Naive_SE_tmp.df; 
    GLM.X  = X0;
    GLM.C  = glmcont;
    GLM.EDtype = EDtype;
    GLM.EDFreq = BCl; 
    
    SPEC.X_RES   = resNaiveSXp;
    SPEC.Y_RES   = resNaiveYp;
    SPEC.X_pwRES = dpwRESXp;
    SPEC.Y_pwRES = dpwRESYp;
    
    PW.dt     = TempTreMethod;
    PW.dtl    = NumTmpTrend;
    PW.pwmeth = pwdmethod;
    PW.fwhm   = lFWHM;
    PW.MAp    = MPparamNum;
    PW.ARp    = Mord;
    PW.nonSPD = nonstationaryvox;
    
    MatFileName = [Path2ImgResults '/sub-' SubID '_ses-' SesID '_ED' EDtype '_' num2str(BCl) '_' pwdmethod '_AR' num2str(Mord) '_MA' num2str(MPparamNum) '_FWHM' num2str(lFWHM) '_' TempTreMethod num2str(NumTmpTrend) '.mat'];
    save(MatFileName,'GLM','SPEC','PW')
end


