clear

COHORT          = 'NEO'; 
lFWHM           = 0; 
icaclean        = 0; 
EDtype          = 'ERF'; 
TempTreMethod   = 'hpf';

MParamNum           = 0;
gsrflag             = 0;
TempDerv            = 0;

disp('=======================================')

PATH2AUX='~/bin/FILM2';
addpath([PATH2AUX '/utils/Trend'])
addpath([PATH2AUX '/utils/AR_YW'])
addpath([PATH2AUX '/mis'])
addpath (fullfile ('/users/nichols/scf915', 'spm12-r7771'));
%addpath('/well/nichols/users/scf915/externals/spm12')


disp('=====SET UP PATHS =============================')

%Path2ImgRaw=[PATH2AUX '/ExampleData/R.mpp'];
%Path2ImgDir = ['/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/ROCKLAND/sub-' SubID '/ses-' SesID '/sub-' SubID '_ses-' SesID '_task-rest_acq-645_bold_mpp'];
%Path2ImgDir = ['/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/' COHORT '/' SubID '_3T_rfMRI_' SesID '_mpp/'];
COHORTDIR='/well/nichols/users/kfh142/data/baby/neofmri_2nd_release_rerun2/';

SubList = {'CC00649XX23','CC00698XX23','CC00789XX23','CC00797XX23','CC00839XX23','CC00847XX23'};
SesList = {'191201','220400','21110','12110','23710','26910'};


for s_cnt = 1:numel(SubList)
    
    SubID = SubList{s_cnt};
    SesID = SesList{s_cnt};
    
    % What is flowing in from the cluster:
    disp('From the cluster ======================')
    disp(['SubID: ' SubID])
    disp(['SesID: ' SesID])
    
    if strcmpi(COHORT,'ROCKLAND')
        Path2ImgRaw = [COHORTDIR '/R_mpp/sub-' SubID '/ses-' SesID];
        Path2ImgDir = [Path2ImgRaw '/sub-' SubID '_ses-' SesID '_task-rest_acq-' num2str(TR*1000) '_bold_mpp'];
        SEGmaskinEPI = [Path2ImgDir '/seg/sub-' SubID '_ses-' SesID '_T1w_func_seg.nii.gz'];
        WMseg        = [Path2ImgDir '/seg/sub-' SubID '_ses-' SesID '_T1w_func_seg_2.nii.gz'];
        path2mask    = [Path2ImgDir '/mask.nii.gz']; 
        Path2MC      = [Path2ImgDir '/prefiltered_func_data_mcf.par'];

    elseif strcmpi(COHORT,'HCP')
        Path2ImgRaw = [COHORTDIR '/R_mpp/sub-' SubID '/ses-' SesID];
        Path2ImgDir = [Path2ImgRaw '/' SubID '_3T_rfMRI_' SesID '_mpp'];
        SEGmaskinEPI = [Path2ImgDir '/seg/' SubID '_3T_T1w_MPR1_func_seg.nii.gz'];
        WMseg        = [Path2ImgDir '/seg/' SubID '_3T_T1w_MPR1_func_seg_2.nii.gz'];
        path2mask    = [Path2ImgDir '/mask.nii.gz']; 
        Path2MC      = [Path2ImgDir '/prefiltered_func_data_mcf.par'];

    elseif any(strcmpi(COHORT,{'Beijing','Cambridge'}))
        Path2ImgRaw  = [COHORTDIR '/R_mpp/sub-' SubID '/ses-' SesID];
        Path2ImgDir  = [Path2ImgRaw '/rest_mpp'];
        SEGmaskinEPI = [Path2ImgDir '/seg/mprage_T1_func_seg.nii.gz'];
        WMseg        = [Path2ImgDir '/seg/mprage_T1_func_seg_2.nii.gz'];
        path2mask    = [Path2ImgDir '/mask.nii.gz']; 
        Path2MC      = [Path2ImgDir '/prefiltered_func_data_mcf.par'];

    elseif strcmpi(COHORT,'NEO')
        Path2ImgRaw  = [COHORTDIR   '/sub-' SubID '/ses-' SesID];
        Path2ImgDir  = Path2ImgRaw;
        SEGmaskinEPI = [Path2ImgDir '/fix/func_dseg_seg.nii.gz'];
        WMseg        = [Path2ImgDir '/fix/func_dseg_wm.nii.gz']; 
        path2mask    = [Path2ImgDir '/fix/func_dseg_mask.nii.gz']; 
        Path2MC      = [Path2ImgDir '/mcdc/func_mcdc.eddy_parameters'];    
    end


    fwhmlab='';
    if lFWHM
        fwhmlab=['_fwhm' num2str(lFWHM)];
    end

    if any(strcmpi(COHORT,{'Beijing','Cambridge','HCP','ROCKLAND'}))
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
    elseif strcmpi(COHORT,'NEO')
        if ~icaclean
            icalab = 'off';
            Path2Img    = [Path2ImgDir '/mcdc/func_mcdc_masked_brain' fwhmlab '.nii.gz'];
        elseif icaclean==1
            icalab = 'nonaggr';
            Path2Img    = [Path2ImgDir '/fix/func_fix_masked_brain' fwhmlab '.nii.gz'];
        end

    else
        disp('XXXX UNRECOGNISED DATASET XXXX')
    end

    disp(['Image: ' Path2Img])
    disp(['Motion params: ' Path2MC])

    % Directory 2 save the results
    Path2ImgResults='/well/nichols/users/scf915/NEO/R.PW/ImgSpecCheck';
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
        BCl = 10; 
        [EDX1,EDX2] = Generate2ER(T,TR,1,BCl);
        EDX = [EDX1,EDX2]; 
    elseif strcmpi(EDtype,'erf')
        BCl = 0;
        path2evs=[PATH2AUX '/mis/EVs/' COHORT '/' COHORT '_ERF_T' num2str(T) '_TR' num2str(TR*1000) '.txt'];
        EDX = load(path2evs);    
        disp(['Paradigm comes from: ' path2evs])
    elseif strcmpi(EDtype,'RegpEV')
        BCl = 0;
        path2evs=[PATH2AUX '/mis/EVs/RegpEV_' COHORT '/RegpEV_' COHORT '_sub_' SubID '_T' num2str(T) '_TR' num2str(TR*1000) '.txt'];
        EDX = load(path2evs);    
        disp(['Paradigm comes from: ' path2evs])
        nEDX = size(EDX,2);
    elseif strcmpi(EDtype,'ORPE')
        BCl  = 0;
        EDX  = GenerateORPE(T,TR,1);
        nEDX = size(EDX,2);
        disp('Paradigm was set to be One Regressor Per Event.')
    end
    disp(['The paradigm is: ' EDtype])

    X   = EDX;

    disp(['design updated, ' num2str(size(X,2))])

    if TempDerv
        %add the temporal derivatives
        tdEDX = EDX;
        tdEDX = tdEDX(1:end-1,:)-tdEDX(2:end,:);
        tdEDx = [tdEDX(1,:); tdEDX];
        % 
        X   = [X,tdEDx];
        disp(['design updated, ' num2str(size(X,2))])

    end

    % Motion parameters ----------------------------------------
    if icaclean==-1; MParamNum = 0; end; 

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
    
    
    [Xs,Ys]        = DrawMeSpectrum(Y,TR,0);
    YSpecs(s_cnt,:)   = mean(Ys,2);    
    
    % Temporal trends ----------------------------------------
    TempTrend = []; K = []; 
    if ~exist('NumTmpTrend','var'); NumTmpTrend=[]; end;
    if strcmpi(TempTreMethod,{'hpf'})
        if isempty(NumTmpTrend) || ~exist('NumTmpTrend','var'); NumTmpTrend=100; end; 
        hp_ff = hp_fsl(T,NumTmpTrend,TR);    
        dY     = hp_ff*Y;  % high pass filter the data
    elseif strcmpi(TempTreMethod,{'hpfc'})    
        NumTmpTrend = hpf_cutoffcalc(X,TR,[Path2ImgResults '/' EDtype '_' num2str(BCl) 'design.mat']);
        hp_ff = hp_fsl(T,NumTmpTrend,TR);    
        dY     = hp_ff*Y;  % high pass filter the data      
    elseif strcmpi(TempTreMethod,{'hpfk'})    
        NumTmpTrend = hpf_cutoffcalc(X,TR,[Path2ImgResults '/' EDtype '_' num2str(BCl) 'design.mat']);
        hp_ff = hp_fsl(T,NumTmpTrend,TR);    
        dY     = hp_ff*Y;  % high pass filter the data      
    end
    clear Y 
    

    [dXs,dYs]        = DrawMeSpectrum(dY,TR,0);
    YdSpecs(s_cnt,:)   = mean(dYs,2);
    
    % Global Signal -----------------------------------------
    [dXs,dgYs]        = DrawMeSpectrum(mean(dY,2),TR,0);
    gYdSpecs(s_cnt,:) = dgYs;    
    
    % Fit a dumb GLM ----------------------------------------
    X              = [ones(T,1),X];
    glmcont        = zeros(1,size(X,2));    
    glmcont([2,3]) = [1 -1];
    [~,~,RES,~]    = myOLS(dY,X,glmcont);
    
    [dXs,dYres]         = DrawMeSpectrum(RES,TR,0);
    RESdSpecs(s_cnt,:)  = mean(dYres,2);
    
    clear dY RES X
end

MatFileName = [Path2ImgResults '/ED' EDtype '_' num2str(BCl) '_FWHM' num2str(lFWHM) '_' TempTreMethod num2str(NumTmpTrend) '_MP' num2str(MParamNum) '_ICACLEAN' num2str(icaclean) '_GSR' num2str(gsrflag) '.mat'];
save(MatFileName,'dXs','YSpecs','YdSpecs','gYdSpecs','RESdSpecs')

disp('xxDONExx')
