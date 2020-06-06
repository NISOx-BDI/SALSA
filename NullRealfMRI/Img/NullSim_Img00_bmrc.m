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
disp(['COHORT name:' COHORT])
disp(['COHORT directory:' COHORTDIR])
disp(['Parth 2 Results: ' Path2ImgResults ])
disp(['FSLDIR: ' getenv('FSLDIR')])
disp(['EDtype: ' EDtype])
disp('=======================================')


SaveImagesFlag      = 1; 
SaveMatFileFlag     = 1; 
MParamNum           = 24;
%gsrflag             = 1;
%icaclean            = 2;

PATH2AUX='~/bin/FILM2';
addpath([PATH2AUX '/utils/Trend'])
addpath([PATH2AUX '/utils/AR_YW'])
addpath([PATH2AUX '/utils/ARMA_HR'])
addpath([PATH2AUX '/mis'])
addpath (fullfile ('/users/nichols/scf915', 'spm12-r7771'));
%addpath('/well/nichols/users/scf915/externals/spm12')

disp('=====SET UP PATHS =============================')
%Raw Images (MMP feat output)
Path2ImgRaw = [COHORTDIR '/R_mpp/sub-' SubID '/ses-' SesID];
if strcmpi(COHORT,'ROCKLAND')
    Path2ImgDir = [Path2ImgRaw '/sub-' SubID '_ses-' SesID '_task-rest_acq-' num2str(TR*1000) '_bold_mpp'];
elseif any(strcmpi(COHORT,{'Beijing','Cambridge'}))
    Path2ImgDir = [Path2ImgRaw '/rest_mpp'];
end


fwhmlab='';
if lFWHM
    fwhmlab=['_fwhm' num2str(lFWHM)];
end

if ~icaclean
    icalab = 'off';
    Path2Img    = [Path2ImgDir '/prefiltered_func_data_bet' fwhmlab '.nii.gz'];
elseif icaclean==1
    icalab = 'nonaggr';
    Path2Img    = [Path2ImgDir '/ica-aroma' fwhmlab '/denoised_func_data_nonaggr.nii.gz'];
elseif icaclean==2
    icalab = 'aggr';
    Path2Img    = [Path2ImgDir '/ica-aroma' fwhmlab '/denoised_func_data_nonaggr.nii.gz'];
end

path2mask = [Path2ImgDir '/mask.nii.gz']; 
Path2MC   = [Path2ImgDir '/prefiltered_func_data_mcf.par'];

SEGmaskinEPI = [Path2ImgDir '/seg/sub-' SubID '_ses-' SesID '_T1w_func_seg.nii.gz'];
WMseg        = [Path2ImgDir '/seg/sub-' SubID '_ses-' SesID '_T1w_func_seg_2.nii.gz'];

disp(['Image: ' Path2Img])
disp(['Motion params: ' Path2MC])

% Directory 2 save the results
Path2ImgResults=[Path2ImgResults '/' SubID '_' SesID];
if ~exist(Path2ImgResults, 'dir')
	mkdir(Path2ImgResults)
	disp(['The directory: ' Path2ImgResults ' did not exists. I made one. '])
end

disp(['Output stuff: ' Path2ImgResults])

%%% Read The Data %%%%%%%%%%%%%%%%%%%%%%%%
disp('=====LOAD THE IMAGE ===========================')

[Y,InputImgStat]=CleanNIFTI_spm(Path2Img,'demean','norm',100);

%---------------------------------------------

T     = InputImgStat.CleanedDim(2);
TR    = InputImgStat.voxelsize(4);
Vorig = InputImgStat.CleanedDim(1);
V     = Vorig;

if size(Y,1)~=T; Y = Y'; end; %TxV

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DESIGN MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('++++++++++++++++++++++++++++++++++++')
disp('++++++++++++ Construct a design matrix')
disp('++++++++++++++++++++++++++++++++++++')

%%% Generate a Design Matrix --------------------------------

if strcmpi(EDtype,'boxcar')
    BCl = 20;
    EDX = GenerateED(BCl,T,TR,fix(T/15)); 
    EDX = EDX - mean(EDX); 
elseif strcmpi(EDtype,'er')
    BCl = 0;
    path2evs=[PATH2AUX '/mis/EVs/' COHORT '/' COHORT '_sub_' SubID '_T' num2str(T) '_TR' num2str(TR*1000) '.txt'];
    EDX = load(path2evs);
    disp(['Paradigm comes from: ' path2evs])
elseif strcmpi(EDtype,'erf')
    BCl = 0;
    path2evs=[PATH2AUX '/mis/EVs/' COHORT '/' COHORT '_ERF_T' num2str(T) '_TR' num2str(TR*1000) '.txt'];
    EDX = load(path2evs);    
    disp(['Paradigm comes from: ' path2evs])
end
disp(['The paradigm is: ' EDtype])

X   = EDX;

disp(['design updated, ' num2str(size(X,2))])

%add the temporal derivatives
tdEDX = EDX;
tdEDX = tdEDX(1:end-1,:)-tdEDX(2:end,:);
tdEDx = [tdEDX(1,:); tdEDX];

X   = [EDX tdEDx];
disp(['design updated, ' num2str(size(X,2))])
% Motion parameters ----------------------------------------
MCp = load(Path2MC);
MCp = GenMotionParam(MCp,MParamNum); 

X = [X,MCp];
disp(['design updated, ' num2str(size(X,2))])

% Global Signal -----------------------------------------
if gsrflag
    GSRts = mean(Y,2); 
    X = [X,GSRts];
    disp(['global signal regression: ' num2str(size(GSRts,1)) ' x ' num2str(size(GSRts,2))])
    disp(['design updated, ' num2str(size(X,2))]) 
end

% Temporal trends ----------------------------------------
TempTrend = []; K = []; 
if ~exist('NumTmpTrend','var'); NumTmpTrend=[]; end;
if any(strcmpi(TempTreMethod,{'dct','spline','poly'}))
    [TempTrend,NumTmpTrend]   = GenerateTemporalTrends(T,TR,TempTreMethod,NumTmpTrend); % DC + Poly trends + Spline trends 
    TempTrend   = TempTrend(:,2:end); % we add a column of one later.
elseif strcmpi(TempTreMethod,{'hpf'})
    if isempty(NumTmpTrend) || ~exist('NumTmpTrend','var'); NumTmpTrend=100; end; 
    hp_ff = hp_fsl(T,NumTmpTrend,TR);    
    X     = hp_ff*X;    % high pass filter the design
    Y     = hp_ff*Y;  % high pass filter the data
elseif strcmpi(TempTreMethod,{'hpfc'})    
    NumTmpTrend = hpf_cutoffcalc(X,TR,[Path2ImgResults '/' EDtype '_' num2str(BCl) 'design.mat']);
    hp_ff = hp_fsl(T,NumTmpTrend,TR);    
    X     = hp_ff*X;    % high pass filter the design
    Y     = hp_ff*Y;  % high pass filter the data      
elseif strcmpi(TempTreMethod,{'hpfk'})    
    NumTmpTrend = hpf_cutoffcalc(X,TR,[Path2ImgResults '/' EDtype '_' num2str(BCl) 'design.mat']);
    hp_ff = hp_fsl(T,NumTmpTrend,TR);    
    X     = hp_ff*X;    % high pass filter the design
    Y     = hp_ff*Y;  % high pass filter the data      
    K     = hp_ff;
end
disp(['Detrending: ' TempTreMethod ',param: ' num2str(NumTmpTrend)])
%
X           = [X,TempTrend];
disp(['design updated, ' num2str(size(X,2))])

% Centre the design  ----------------------------------
X           = X - mean(X); % demean everything 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIT A MODEL TO THE ACd DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X                  = [ones(T,1),X];
glmcont            = zeros(1,size(X,2));

if strcmpi(EDtype,'boxcar')
    glmcont(2)          = 1; 
    disp('+ single contrast for boxcar is set.')
elseif strcmpi(EDtype,'er') || strcmpi(EDtype,'erf')
    glmcont([2,3])     = [1 -1];
    disp('+ double tasks contrasted against each other.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PREWHITEN THE RESIDULAS & ESTIMATE BIAS AND CPS %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('++++++++++++++++++++++++++++++++++++')
disp('++++++++++++PREWHITEN THE MODEL.')
disp('++++++++++++++++++++++++++++++++++++')
MPparamNum = 0; 
if strcmpi(pwdmethod,'AR-W') %Worsely %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,WBLUSRES,wse,wtv,wzv] = arw(Y,X,glmcont,Mord,InputImgStat,path2mask,K);
elseif strcmpi(pwdmethod,'gFAST') %SPMfast %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,WBLUSRES,wse,wtv,wzv] = gfast(Y,X,TR,glmcont,1);
elseif strcmpi(pwdmethod,'ACFadj') % Yule-Walker %%%%%%%%%%%%%%%%%%%%%%%%%%
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,WBLUSRES,wse,wtv,wzv] = feat5(Y,X,glmcont,Mord,InputImgStat,path2mask,1,0,K);   
elseif strcmpi(pwdmethod,'gFASTxACFadj') % Yule-Walker %%%%%%%%%%%%%%%%%%%%    
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,WBLUSRES,wse,wtv,wzv] = gReMLxACF(Y,X,TR,glmcont,Mord,InputImgStat,path2mask,1,K);
elseif strcmpi(pwdmethod,'gFASTxACFadj2t') % Yule-Walker %%%%%%%%%%%%%%%%%%%%    
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,WBLUSRES,wse,wtv,wzv] = gReMLxACF(Y,X,TR,glmcont,Mord,InputImgStat,path2mask,1,K,WMseg);    
elseif strcmpi(pwdmethod,'gFASTxACFadj2t2j') % Yule-Walker %%%%%%%%%%%%%%%%%%%%    
    J = 2; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,WBLUSRES,wse,wtv,wzv] = gReMLxACF(Y,X,TR,glmcont,Mord,InputImgStat,path2mask,1,K,WMseg,J);        
elseif strcmpi(pwdmethod,'ACF') % ACF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,WBLUSRES,wse,wtv,wzv] = feat5(Y,X,glmcont,Mord,InputImgStat,path2mask,0,[]);
elseif strcmpi(pwdmethod,'ARMAHR') % ARMAHR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error('Use the old script.')
elseif strcmpi(pwdmethod,'ARMAReML') % ARMAHR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error('Use the old script.')
else
    error('Unrecognised prewhitening method.')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACL OF THE RESIDUALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('++++++++++++Calculate the ACL')
acl  = ACLImage(RES');
wacl = ACLImage(WRES');

disp(['mean acl of res   : ' num2str(mean(acl)) ])
disp(['mean acl of w res : ' num2str(mean(wacl)) ])

disp('+++++++++++Calculate AR1')
ar1  = AR_YW_voxel(RES,T,1);
war1 = AR_YW_voxel(WRES,T,1);

disp(['mean ar1 of res   : ' num2str(mean(ar1)) ])
disp(['mean ar1 of w res : ' num2str(mean(war1)) ])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('++++++++++++++++++++++++++++++++++++')
disp('++++++++++++Calculate the CPS.')
disp('++++++++++++++++++++++++++++++++++++')

% CPS on naive
%BLUSRES   = BLUSres(Y,X,1:size(X,2)); % temp off, too much exec time. 
[~,~,cpz] = CPSUnivar(RES,stat.df);

%CPS on whitened residuals
[~,~,wcpz] = CPSUnivar(WRES,stat.df);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECTRUM OF THE RESIDUALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('++++++++++++++++++++++++++++++++++++')
disp('++++++++++++Calculate the spectrum of the residuals.')
disp('++++++++++++++++++++++++++++++++++++')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Whole brain
[dpwRESXp,dpwRESYp]             = DrawMeSpectrum(WRES,TR,0);
dpwRESYp                        = mean(dpwRESYp,2); % average across voxels

% Mask the residuals
WRES                            = MaskImg(WRES',SEGmaskinEPI,InputImgStat);

% Grey Matter
[dpwRESXp_GM,dpwRESYp_GM]       = DrawMeSpectrum(WRES{2}',TR,0);
dpwRESYp_GM                     = mean(dpwRESYp_GM,2); % average across voxels

% White Matter
[dpwRESXp_WM,dpwRESYp_WM]       = DrawMeSpectrum(WRES{3}',TR,0);
dpwRESYp_WM                     = mean(dpwRESYp_WM,2); % average across voxels

%CSF
[dpwRESXp_CSF,dpwRESYp_CSF]     = DrawMeSpectrum(WRES{1}',TR,0);
dpwRESYp_CSF                    = mean(dpwRESYp_CSF,2); % average across voxels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear dpwRES
[resNaiveXp,resNaiveYp]        = DrawMeSpectrum(RES,TR,0);
resNaiveYp                      = mean(resNaiveYp,2); % average across voxels

% Mask the residuals
RES                             = MaskImg(RES',SEGmaskinEPI,InputImgStat);

% Grey Matter
[resNaiveXp_GM,resNaiveYp_GM]   = DrawMeSpectrum(RES{2}',TR,0);
resNaiveYp_GM                   = mean(resNaiveYp_GM,2); % average across voxels

% White Matter
[resNaiveXp_WM,resNaiveYp_WM]   = DrawMeSpectrum(RES{3}',TR,0);
resNaiveYp_WM                   = mean(resNaiveYp_WM,2); % average across voxels

%CSF
[resNaiveXp_CSF,resNaiveYp_CSF] = DrawMeSpectrum(RES{1}',TR,0);
resNaiveYp_CSF                  = mean(resNaiveYp_CSF,2); % average across voxels

%clear RES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE THE RESULTS AS AN IMAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('++++++++++++++++++++++++++++++++++++')
disp('++++++++++++Save the results.')
disp('++++++++++++++++++++++++++++++++++++')

if SaveImagesFlag
    % 3D IMAGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    VariableList = {'cbhat','se','tv',...
        'Wcbhat','wse','wtv',...
        'acl','wacl',...
        'cpz','wcpz',...
        'ar1','war1'};
    OutputImgStat            = InputImgStat.spmV(1);
    OutputImgStat.Removables = InputImgStat.Removables;

    for vname = VariableList

        tmpvar     = eval(vname{1});
        fname      = [Path2ImgResults '/ED' EDtype '_' num2str(BCl) '_' pwdmethod '_AR' num2str(Mord) '_MA' num2str(MPparamNum) '_FWHM' num2str(lFWHM) '_' TempTreMethod num2str(NumTmpTrend) '_' vname{1} '_ICACLEAN' num2str(icaclean) '_GSR' num2str(gsrflag) '.nii'];

        CleanNIFTI_spm(tmpvar,'ImgInfo',InputImgStat.spmV,'DestDir',fname,'removables',InputImgStat.Removables);
        %system(['gzip ' fname]);
    end
end

% MAT FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SaveMatFileFlag
    GLM.df      = stat.df; 
    GLM.X       = X;
    GLM.C       = glmcont;
    GLM.EDtype  = EDtype;
    GLM.EDFreq  = BCl; 
    
    % Overal Specturm
    SPEC.X_RES   = resNaiveXp;
    SPEC.Y_RES   = resNaiveYp;
    SPEC.X_pwRES = dpwRESXp;
    SPEC.Y_pwRES = dpwRESYp;
    
    % Spectrum of GM
    SPEC.X_RES_GM   = resNaiveXp_GM;
    SPEC.Y_RES_GM   = resNaiveYp_GM;
    SPEC.X_pwRES_GM = dpwRESXp_GM;
    SPEC.Y_pwRES_GM = dpwRESYp_GM;
    
    % Spectrum of WM
    SPEC.X_RES_WM   = resNaiveXp_WM;
    SPEC.Y_RES_WM   = resNaiveYp_WM;
    SPEC.X_pwRES_WM = dpwRESXp_WM;
    SPEC.Y_pwRES_WM = dpwRESYp_WM;
    
    % Spectrum of CSF
    SPEC.X_RES_CSF   = resNaiveXp_CSF;
    SPEC.Y_RES_CSF   = resNaiveYp_CSF;
    SPEC.X_pwRES_CSF = dpwRESXp_CSF;
    SPEC.Y_pwRES_CSF = dpwRESYp_CSF;    
    
    PW.dt     = TempTreMethod;
    PW.dtl    = NumTmpTrend;
    PW.pwmeth = pwdmethod;
    PW.fwhm   = lFWHM;
    PW.MAp    = MPparamNum;
    PW.ARp    = Mord;
    
    MatFileName = [Path2ImgResults '/ED' EDtype '_' num2str(BCl) '_' pwdmethod '_AR' num2str(Mord) '_MA' num2str(MPparamNum) '_FWHM' num2str(lFWHM) '_' TempTreMethod num2str(NumTmpTrend) '_ICACLEAN' num2str(icaclean) '_GSR' num2str(gsrflag) '.mat'];
    save(MatFileName,'GLM','SPEC','PW')
end

disp('xxDONExx')