clear; 
%warning('on','all')

pwdmethod = 'ACF'; %ACF AR-YW AR-W ARMAHR
Mord      = 20; 
lFWHM     = 0;
SubID     = 'A00028858';
SesID     = 'DS2'; 
TR        = 0.645;

%TempTreMethod = 'spline'; 
%NumTmpTrend   = 3;

TempTreMethod = 'hpf'; 
NumTmpTrend   = [];





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
%Path2ImgDir = ['/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/R_test/sub-' SubID '_ses-' SesID '_task-rest_acq-645_bold_mpp'];
Path2ImgDir = ['/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/ROCKLAND/sub-' SubID '/ses-' SesID '/sub-' SubID '_ses-' SesID '_task-rest_acq-645_bold_mpp'];
if ~lFWHM
    Path2Img    = [Path2ImgDir '/prefiltered_func_data_gm_bet.nii.gz'];
else
    Path2Img    = [Path2ImgDir '/prefiltered_func_data_bet_FWHM' num2str(lFWHM) '.nii'];
end

Path2MC  = [Path2ImgDir '/prefiltered_func_data_mcf.par'];

disp(['Image: ' Path2Img])
disp(['Motion params: ' Path2MC])

% Directory 2 save the results
Path2ImgResults=[PATH2AUX '/ExampleData/R.mpp/RNullfMRI_' SubID '_' SesID];
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
% ----------------------------------------------

CLK	 = fix(clock);
tmpdir  = [tempdir 'octspm12/tmp_' num2str(randi(5000)) '_' num2str(CLK(end))]; % make a temp directory 
mkdir(tmpdir)
disp(['Unzip into: ' tmpdir ])
randtempfilename=[tmpdir '/prefilt_tmp_' SubID '_' num2str(randi(50)) num2str(CLK(end)+randi(10)) '.nii'];
system(['gunzip -c ' Path2Img ' > ' randtempfilename]); %gunzip function in Octave deletes the source file in my version!

[Y,InputImgStat]=CleanNIFTI_spm(randtempfilename,'demean');

disp(['Remove the temp directory: ' tmpdir])
%rmdir(tmpdir,'s')
system(['rm -rf ' tmpdir])
%-----------------------------------------------

T = InputImgStat.CleanedDim(2);
TR = InputImgStat.voxelsize(4);
Vorig = InputImgStat.CleanedDim(1);
V = Vorig;


if size(Y,1)~=T; Y = Y'; end; %TxV

Y = Y(:,1:25:end);
[~,V] = size(Y); 
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
    dY    = dY - mean(dY); 
end
disp(['Detrending: ' TempTreMethod ',param: ' num2str(NumTmpTrend)])
%
X           = [X,TempTrend];
disp(['design updated, ' num2str(size(X,2))])

% Centre the design  ----------------------------------
X           = X - mean(X); % demean everything 

X0          = [ones(T,1),X];

%%% RESIDUALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('++++++++++++Get the residuals using I-XX+')
pinvX           = pinv(X0); 
ResidFormingMat = eye(T)-X0*pinvX; % residual forming matrix 
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
[Bhat_Naive,~,resNaive,Stat_Naive_SE_tmp]   = myOLS(dY,X0,glmcont);
SE_Naive                                    = Stat_Naive_SE_tmp.se;
tVALUE_Naive                                = Stat_Naive_SE_tmp.tval;
[~,CPSstat_Naive,CPZ_Naive]                 = CPSUnivar(resNaive,X0);


pw_residY   = PrewhitenSPD(residY,15,1/TR); % whiten the residuals
pw_dY       = Prewhiten2x2SPD(dY,residY,15,1/TR); % whiten the Y with residuals


% for ii = 1:size(X,2)
%     ii
%     pw_X_tmp(ii,:,:) = Prewhiten2x2SPD(X(:,ii),residY,15,1/TR); % whiten the X with residuals
% end

for iv = 1:V
    if ~mod(iv,100); disp(iv); end; 
    Xnew = Prewhiten2x2SPD(X,residY(:,iv),15,1/TR);
    [~,~,dpwRES(:,iv)] = myOLS(pw_dY(:,iv),Xnew);
    
end


[RESxx,RESyy]         = DrawMeSpectrum(dpwRES,TR);
plot(RESxx,mean(RESyy,2))

[dYxx,dYyy]         = DrawMeSpectrum(dY,TR);
[xx,yy]             = DrawMeSpectrum(resNaive,TR);
[pwxx,pwyy]         = DrawMeSpectrum(pw_residY,TR);
[pwYxx,pwYyy]       = DrawMeSpectrum(pw_dY,TR);
%[pwYxx_tt,pwYyy_tt] = DrawMeSpectrum(pw_residY_test,TR);


figure; hold on; grid on; 
plot(dYxx,mean(dYyy,2))
plot(pwxx,mean(pwyy,2))
plot(pwYxx,mean(pwYyy,2))
%plot(pwYxx_tt,mean(pwYyy_tt,2))
plot(xx,mean(yy,2))


