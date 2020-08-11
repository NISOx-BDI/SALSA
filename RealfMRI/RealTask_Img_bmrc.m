
warning('off','all')


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COHORT = 'ROCKLAND'; 
% TaskName = 'CHECKERBOARD'; 
% pwdmethod = 'gACFadjxACFadj2tT1S0P5'; %ACF AR-YW AR-W ARMAHR
% Mord      = 16; 
% lFWHM     = 5;
% SubID     = 'A00008399';
% SesID     = 'DS2'; 
% TR        = 0.645;
% icaclean  = 0;
% %TempTreMethod = 'spline'; 
% %NumTmpTrend   = 3;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/R_ROCKLAND/R.PW/T900/ROCKLAND_645_900_ACFadjT0S0_AR--2_MA-0_FWHM5_hpfk_gsr0_aroma0/A00028150_DS

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
disp('=======================================')




tic;
SaveImagesFlag      = 1; 
SaveMatFileFlag     = 1; 
DoDetrendingPrior   = 0; 
MParamNum           = 24;
gsrflag             = 0;
TempDerv            = 1;

disp('=======================================')

PATH2AUX='~/bin/FILM2';
addpath([PATH2AUX '/utils/Trend'])
addpath([PATH2AUX '/utils/AR_YW'])
addpath([PATH2AUX '/utils/ARMA_HR'])
addpath([PATH2AUX '/mis'])
addpath (fullfile ('/users/nichols/scf915', 'spm12-r7771'))


disp('=====SET UP PATHS =============================')


if strcmpi(COHORT,'ROCKLAND')
    %task-CHECKERBOARD_acq-645
    EDtype = ['task-' TaskName '_acq-' num2str(TR*1000)]; 
    Path2ImgRaw = [COHORTDIR '/R_mpp/sub-' SubID '/ses-' SesID];
    Path2ImgDir = [Path2ImgRaw '/sub-' SubID '_ses-' SesID '_' EDtype '_bold_mpp'];
    SEGmaskinEPI = [Path2ImgDir '/seg/sub-' SubID '_ses-' SesID '_T1w_func_seg.nii.gz'];
    WMseg        = [Path2ImgDir '/seg/sub-' SubID '_ses-' SesID '_T1w_func_seg_2.nii.gz'];
    path2mask    = [Path2ImgDir '/mask.nii.gz']; 
    Path2MC      = [Path2ImgDir '/prefiltered_func_data_mcf.par'];
    
    Path2Event = [Path2ImgDir '/sub-' SubID '_ses-' SesID '_' EDtype '_events.tsv']; 

elseif strcmpi(COHORT,'tHCP')
    TaskName = SesID; 
    EDtype = ['task-' TaskName '_acq-' num2str(TR*1000)];
    Path2ImgRaw = [COHORTDIR '/R_mpp/sub-' SubID '/ses-' SesID];
    Path2ImgDir = [Path2ImgRaw '/' SubID '_3T_tfMRI_' SesID '_mpp'];
    SEGmaskinEPI = [Path2ImgDir '/seg/' SubID '_3T_T1w_MPR1_func_seg.nii.gz'];
    WMseg        = [Path2ImgDir '/seg/' SubID '_3T_T1w_MPR1_func_seg_2.nii.gz'];
    path2mask    = [Path2ImgDir '/mask.nii.gz']; 
    Path2MC      = [Path2ImgDir '/prefiltered_func_data_mcf.par']; 
else
    error('Unrecog data.')
end


fwhmlab='';
if lFWHM
    fwhmlab=['_fwhm' num2str(lFWHM)];
end

if any(strcmpi(COHORT,{'Beijing','Cambridge','tHCP','ROCKLAND'}))
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

[EDX,EventTrail,OnSet,Duration,Events] = readBIDSEvent(Path2Event,T,TR);

X = EDX;

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
    NumTmpTrend = hpf_cutoffcalc(X,TR,[Path2ImgResults '/' EDtype '_design.mat']);
    hp_ff = hp_fsl(T,NumTmpTrend,TR);    
    X     = hp_ff*X;    % high pass filter the design
    Y     = hp_ff*Y;  % high pass filter the data      
elseif strcmpi(TempTreMethod,{'hpfk'})    
    NumTmpTrend = hpf_cutoffcalc(X,TR,[Path2ImgResults '/' EDtype '_design.mat']);
    hp_ff = hp_fsl(T,NumTmpTrend,TR);    
    X     = hp_ff*X;    % high pass filter the design
    Y     = hp_ff*Y;  % high pass filter the data      
    K     = hp_ff;
end
disp(['Detrending: ' TempTreMethod ',param: ' num2str(NumTmpTrend)])
%
X           = [X,TempTrend];
disp(['design updated, ' num2str(size(X,2))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIT A MODEL TO THE ACd DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Centre the design  ----------------------------------
X           = X - mean(X); % demean everything 
X           = [ones(T,1),X];
glmcont     = zeros(1,size(X,2));
if strcmpi(TaskName,'CHECKERBOARD') && strcmpi(COHORT,'ROCKLAND')
    glmcont(2)          = 1; 
    disp(['+ single contrast for boxcar is set. TASK: ' TaskName])
elseif strcmpi(TaskName,'MOTOR_LR') && strcmpi(COHORT,'tHCP')
    glmcont(2)          = 1; 
    disp('+ single contrast for boxcar is set.')    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PREWHITEN THE RESIDULAS & ESTIMATE BIAS AND CPS %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('++++++++++++++++++++++++++++++++++++')
disp('++++++++++++PREWHITEN THE MODEL.')
disp('++++++++++++++++++++++++++++++++++++')
MPparamNum = 0; 
if strcmpi(pwdmethod,'AR-W') %Worsely %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,~,wse,wtv,wzv] = arw(Y,X,glmcont,Mord,InputImgStat,path2mask,K);
    
% -------------------- FAST     
elseif strcmpi(pwdmethod,'gFAST') %SPMfast %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,~,wse,wtv,wzv] = gfast(Y,X,TR,glmcont,1);
elseif strcmpi(pwdmethod,'vFAST100') %SPMfast %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    aclageval = 100; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = vfast(Y,X,TR,glmcont,Mord,aclageval,K);    
elseif strcmpi(pwdmethod,'vFAST50') %SPMfast %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    aclageval = 50; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = vfast(Y,X,TR,glmcont,Mord,aclageval,K);
elseif strcmpi(pwdmethod,'vFAST20') %SPMfast %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    aclageval = 20; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = vfast(Y,X,TR,glmcont,Mord,aclageval,K);      
elseif strcmpi(pwdmethod,'cFAST') %SPMfast %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = cfast5(Y,X,TR,glmcont,InputImgStat,WMseg);
    
% -------------------- FAST FEAT ------------------------------------------
elseif strcmpi(pwdmethod,'FASTFEAT1') %SPMfast %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     aclageval = 1; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = fastfeat(Y,X,TR,glmcont,Mord,ACFRegF,aclageval,[],[],1,K);
elseif strcmpi(pwdmethod,'FASTFEAT5') %SPMfast %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     aclageval = 5; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = fastfeat(Y,X,TR,glmcont,Mord,ACFRegF,aclageval,[],[],1,K);
elseif strcmpi(pwdmethod,'FASTFEAT10') %SPMfast %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     aclageval = 10; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = fastfeat(Y,X,TR,glmcont,Mord,ACFRegF,aclageval,[],[],1,K);    
elseif strcmpi(pwdmethod,'FASTFEAT20') %SPMfast %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     aclageval = 20; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = fastfeat(Y,X,TR,glmcont,Mord,ACFRegF,aclageval,[],[],1,K);
elseif strcmpi(pwdmethod,'FASTFEAT50') %SPMfast %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     aclageval = 50; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = fastfeat(Y,X,TR,glmcont,Mord,ACFRegF,aclageval,[],[],1,K);
elseif strcmpi(pwdmethod,'FASTFEAT100') %SPMfast %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     aclageval = 100; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = fastfeat(Y,X,TR,glmcont,Mord,ACFRegF,aclageval,[],[],1,K);  
    
elseif strcmpi(pwdmethod,'FASTFEAT20T0S0') %SPMfast %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     aclageval = 20; 
     ACFRegF   = 0; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = fastfeat(Y,X,TR,glmcont,Mord,ACFRegF,aclageval,[],[],1,K);
elseif strcmpi(pwdmethod,'FASTFEAT50T0S0') %SPMfast %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     aclageval = 50;
     ACFRegF   = 0; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = fastfeat(Y,X,TR,glmcont,Mord,ACFRegF,aclageval,[],[],1,K);
elseif strcmpi(pwdmethod,'FASTFEAT100T0S0') %SPMfast %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     aclageval = 100;
     ACFRegF   = 0;
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = fastfeat(Y,X,TR,glmcont,Mord,ACFRegF,aclageval,[],[],1,K);      

    
elseif strcmpi(pwdmethod,'FASTFEAT20T0S5') %SPMfast %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     aclageval = 20; 
     ACFRegF   = 0; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = fastfeat(Y,X,TR,glmcont,Mord,ACFRegF,aclageval,InputImgStat,path2mask,1,K);
elseif strcmpi(pwdmethod,'FASTFEAT50T0S5') %SPMfast %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     aclageval = 50;
     ACFRegF   = 0; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = fastfeat(Y,X,TR,glmcont,Mord,ACFRegF,aclageval,InputImgStat,path2mask,1,K);
elseif strcmpi(pwdmethod,'FASTFEAT100T0S5') %SPMfast %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     aclageval = 100;
     ACFRegF   = 0;
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = fastfeat(Y,X,TR,glmcont,Mord,ACFRegF,aclageval,InputImgStat,path2mask,1,K);          
    
    
% ------------------------ ACFadj -----------------------------------------
elseif strcmpi(pwdmethod,'ACF') % ACF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ACFRegF = 1; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = feat5(Y,X,glmcont,Mord,ACFRegF,InputImgStat,path2mask,0,[]);
elseif strcmpi(pwdmethod,'ACFT1S0') % ACF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ACFRegF = 1; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = feat5(Y,X,glmcont,Mord,ACFRegF,[],[],0,[]);    
elseif strcmpi(pwdmethod,'ACFT0S0') % ACF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ACFRegF = 0; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = feat5(Y,X,glmcont,Mord,ACFRegF,[],[],0,[]);    
elseif strcmpi(pwdmethod,'ACFadj') % Yule-Walker %%%%%%%%%%%%%%%%%%%%%%%%%%
    ACFRegF = 1; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = feat5(Y,X,glmcont,Mord,ACFRegF,InputImgStat,path2mask,1,K); 
elseif strcmpi(pwdmethod,'ACFadjT1S0') % Yule-Walker %%%%%%%%%%%%%%%%%%%%%%%%%%
    ACFRegF      = 1; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = feat5(Y,X,glmcont,Mord,ACFRegF,[],[],1,K);     
elseif strcmpi(pwdmethod,'ACFadjT0S0') % Yule-Walker %%%%%%%%%%%%%%%%%%%%%%%%%%
    ACFRegF      = 0; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = feat5(Y,X,glmcont,Mord,ACFRegF,[],[],1,K);   
elseif strcmpi(pwdmethod,'ACFadjT0S5') % Yule-Walker %%%%%%%%%%%%%%%%%%%%%%%%%%
    ACFRegF      = 0; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = feat5(Y,X,glmcont,Mord,ACFRegF,InputImgStat,path2mask,1,K);   

% --------------------------------------------------------------------------------------
% ------------------------ Global FEAT -------------------------------------------------     
    
elseif strcmpi(pwdmethod,'gACFadjT1S0') % Two stage without tissue segmentation 
    aclageval = 0; 
    ACFRegF   = 1;
    poolflag  = 1;
    [~,~,cbhat,RES,stat,se,tv,zv,Wcbhat,WRES,wse,wtv,wzv] = gfeat(Y,X,TR,glmcont,Mord,ACFRegF,aclageval,1,K,poolflag);

elseif strcmpi(pwdmethod,'gsTEST') % Two stage without tissue segmentation 
    aclageval = 0; 
    ACFRegF   = 1;
    poolflag  = 5;
    [~,~,cbhat,RES,stat,se,tv,zv,Wcbhat,WRES,wse,wtv,wzv] = gsfeat(Y,X,TR,glmcont,Mord,ACFRegF,aclageval,1,K,poolflag);    
    
elseif strcmpi(pwdmethod,'gACFadjT1S0P20') % Two stage without tissue segmentation 
    aclageval = 0; 
    ACFRegF   = 1;
    poolflag  = 20;
    [~,~,~,cbhat,RES,stat,se,tv,zv,Wcbhat,WRES,wse,wtv,wzv] = gsfeat(Y,X,TR,glmcont,Mord,ACFRegF,aclageval,1,K,poolflag);   
    
elseif strcmpi(pwdmethod,'gACFadjT1S0P5AR1') % Two stage without tissue segmentation 
    aclageval = 0; 
    ACFRegF   = 1;
    poolflag  = 3;
    [~,~,~,cbhat,RES,stat,se,tv,zv,Wcbhat,WRES,wse,wtv,wzv] = gsfeat0(Y,X,TR,glmcont,Mord,ACFRegF,aclageval,1,K,poolflag);
    
elseif strcmpi(pwdmethod,'gACFadjT1S0P5') % Two stage without tissue segmentation 
    aclageval = 0; 
    ACFRegF   = 1;
    poolflag  = 5;
    [~,~,~,cbhat,RES,stat,se,tv,zv,Wcbhat,WRES,wse,wtv,wzv] = gsfeat(Y,X,TR,glmcont,Mord,ACFRegF,aclageval,1,K,poolflag);        
% --------------------------------------------------------------------------------------
% ------------------------ Two Stage Methods: gfeatxfeat -------------------------------    

elseif strcmpi(pwdmethod,'gACFadjxACFadjT1S0') % Two stage without tissue segmentation 
    ACFRegF   = 1;
    poolflag  = [];
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = gfeatxfeat(Y,X,TR,glmcont,Mord,ACFRegF,[],[],1,K,[],poolflag);  

elseif strcmpi(pwdmethod,'gACFadjxACFadj2tT1S0') % Two stage with knowledge of tissue
    ACFRegF   = 1; 
    poolflag  = 1; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = gfeatxfeat(Y,X,TR,glmcont,Mord,ACFRegF,InputImgStat,[],1,K,WMseg,poolflag);   
elseif strcmpi(pwdmethod,'gACFadjxACFadj2tT1S0P0') % Two stage with knowledge of tissue
    ACFRegF   = 1; 
    poolflag  = 0; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = gfeatxfeat(Y,X,TR,glmcont,Mord,ACFRegF,InputImgStat,[],1,K,WMseg,poolflag);  
    
elseif strcmpi(pwdmethod,'gACFadjxACFadj2tT1S0P5') % Two stage with knowledge of tissue
    ACFRegF   = 1; 
    poolflag  = 5; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = gsfeatxfeat(Y,X,TR,glmcont,Mord,ACFRegF,InputImgStat,[],1,K,WMseg,poolflag);   
elseif strcmpi(pwdmethod,'gACFadjxACFadjT1S0P5') % Two stage with knowledge of tissue
    ACFRegF   = 1; 
    poolflag  = 5; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = gsfeatxfeat(Y,X,TR,glmcont,Mord,ACFRegF,InputImgStat,[],1,K,[],poolflag);      

% ------------------------ Two stage Methods: gReML x ACF ------------------------------    
elseif strcmpi(pwdmethod,'gFASTxACFadj') % Yule-Walker %%%%%%%%%%%%%%%%%%%%    
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = gReMLxACF(Y,X,TR,glmcont,Mord,InputImgStat,path2mask,1,K);
elseif strcmpi(pwdmethod,'gFASTxACFadj2t') % Yule-Walker %%%%%%%%%%%%%%%%%%%%    
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = gReMLxACF(Y,X,TR,glmcont,Mord,InputImgStat,path2mask,1,K,WMseg);  
elseif strcmpi(pwdmethod,'gFASTxACFadj2t2j') % Yule-Walker %%%%%%%%%%%%%%%%%%%%    
    J            = 2; 
    ACFRegF      = 1; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = gReMLxACF(Y,X,TR,glmcont,Mord,ACFRegF,InputImgStat,path2mask,1,K,WMseg,J);      
elseif strcmpi(pwdmethod,'gFASTxACFadj') % Yule-Walker %%%%%%%%%%%%%%%%%%%%    
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = gReMLxACF(Y,X,TR,glmcont,Mord,InputImgStat,path2mask,1,K);
elseif strcmpi(pwdmethod,'gFASTxACFadj2t') % Yule-Walker %%%%%%%%%%%%%%%%%%%%    
    J            = 2; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = gReMLxACF(Y,X,TR,glmcont,Mord,InputImgStat,path2mask,1,K,WMseg,J);        
elseif strcmpi(pwdmethod,'gFASTxACFadj2t2jT0S0') % Yule-Walker %%%%%%%%%%%%%%%%%%%%    
    J            = 2; 
    path2mask    = [];
    ACFRegF      = 0; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = gReMLxACF(Y,X,TR,glmcont,Mord   ,ACFRegF,InputImgStat,path2mask,1,K,WMseg,J);      
% ------------------------ Two stage Methods: gReML x AR-W ------------------------------        
elseif strcmpi(pwdmethod,'gReMLxARw2t2j') % Yule-Walker %%%%%%%%%%%%%%%%%%%%    
    J = 2; 
    [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = gReMLxARw(Y,X,TR,glmcont,Mord,InputImgStat,path2mask,K,WMseg,J);  

% ------------------------ ARMA MODELS ------------------------------------
elseif strcmpi(pwdmethod,'ARMAHR') % ARMAHR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error('Use the old script.')
elseif strcmpi(pwdmethod,'ARMAReML') % ARMAHR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error('Use the old script.')
else
    error('Unrecognised prewhitening method.')

end
toc
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


wres = WRES; 
res  = RES;
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
    if ~strcmpi(EDtype,'ORPE')
%         VariableList = {'cbhat','se','tv',...
%             'Wcbhat','wse','wtv',...
%             'acl','wacl',...
%             'cpz','wcpz',...
%             'ar1','war1'};
        VariableList = {'tv','wtv'};        
    else
%         VariableList = {'cbhat','se','tv',...
%             'Wcbhat','wse','wtv',...
%             'acl','wacl',...
%             'ar1','war1'};
        VariableList = {'tv','wtv'};        
    end
    
    OutputImgStat            = InputImgStat.spmV(1);
    OutputImgStat.Removables = InputImgStat.Removables;

    for vname = VariableList

        tmpvar     = eval(vname{1});
        fname      = [Path2ImgResults '/ED' EDtype '_' pwdmethod '_AR' num2str(Mord) '_MA' num2str(MPparamNum) '_FWHM' num2str(lFWHM) '_' TempTreMethod num2str(NumTmpTrend) '_' vname{1} '_ICACLEAN' num2str(icaclean) '_GSR' num2str(gsrflag) '.nii'];

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
    
    MatFileName = [Path2ImgResults '/ED' EDtype '_' pwdmethod '_AR' num2str(Mord) '_MA' num2str(MPparamNum) '_FWHM' num2str(lFWHM) '_' TempTreMethod num2str(NumTmpTrend) '_ICACLEAN' num2str(icaclean) '_GSR' num2str(gsrflag) '.mat'];
    save(MatFileName,'GLM','SPEC','PW')
end

disp('xxDONExx')
