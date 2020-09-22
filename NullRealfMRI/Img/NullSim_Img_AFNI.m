clear

addpath('/users/nichols/scf915/bin/FILM2/mis')

COHORT='ROCKLAND';
Tall = 900; 
TR = 0.645;

SubID = 'A00028753';
SesID = 'DS2';

lFWHM = 5; 
BCl = 0; 
EDtype = 'ER';

pwdmethod = '3dREMLfit';
Mord = 1;
MPparamNum = 1;

TempTreMethod = 'poly';
NumTmpTrend = 4;

icaclean = 0;
gsrflag = 0; 

COHORTDIR=['/well/nichols/users/scf915/' COHORT];

load([COHORTDIR '/' COHORT '_sesid.mat'])
load([COHORTDIR '/' COHORT '_subid.mat'])
SubList = participants;
SesList = ses;

for subcnt = 1:51
    T = Tall;
    
    SubID = SubList{subcnt};
    SesID = SesList{subcnt};

    disp([SubID '_' SesID])
    
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

    Path2Dir = [COHORTDIR '/R.PW/' COHORT '_' num2str(TR*1000) '_' num2str(T) '_' pwdmethod '_AR-' num2str(Mord) '_MA-' num2str(MPparamNum) '_FWHM' num2str(lFWHM) '_' TempTreMethod '_gsr' num2str(gsrflag) '_aroma' num2str(icaclean) '/' SubID '_' SesID];
    Path2Img = [Path2Dir '/ER_sub-' SubID '_ses-' SesID '_Decon_Rwherr.nii.gz'];

    if ~exist(Path2Img,'file')
        disp(['DOESNOT EXISTS: ' Path2Img ])
        continue;
    end
    
    [WRES,InputImgStat] = CleanNIFTI_spm(Path2Img);
    T                   = InputImgStat.CleanedDim(2);
    TR                  = InputImgStat.voxelsize(4);


    [dpwRESXp,dpwRESYp]             = DrawMeSpectrum(WRES',TR,0);
    dpwRESYp                        = mean(dpwRESYp,2); % average across voxels

    % Mask the residuals
    WRES                            = MaskImg(WRES,SEGmaskinEPI,InputImgStat);

    % Grey Matter
    [dpwRESXp_GM,dpwRESYp_GM]       = DrawMeSpectrum(WRES{2}',TR,0);
    dpwRESYp_GM                     = mean(dpwRESYp_GM,2); % average across voxels

    % White Matter
    [dpwRESXp_WM,dpwRESYp_WM]       = DrawMeSpectrum(WRES{3}',TR,0);
    dpwRESYp_WM                     = mean(dpwRESYp_WM,2); % average across voxels

    %CSF
    [dpwRESXp_CSF,dpwRESYp_CSF]     = DrawMeSpectrum(WRES{1}',TR,0);
    dpwRESYp_CSF                    = mean(dpwRESYp_CSF,2); % average across voxels


    % Overal Specturm
    SPEC.X_pwRES = dpwRESXp;
    SPEC.Y_pwRES = dpwRESYp;

    % Spectrum of GM
    SPEC.X_pwRES_GM = dpwRESXp_GM;
    SPEC.Y_pwRES_GM = dpwRESYp_GM;

    % Spectrum of WM
    SPEC.X_pwRES_WM = dpwRESXp_WM;
    SPEC.Y_pwRES_WM = dpwRESYp_WM;

    % Spectrum of CSF
    SPEC.X_pwRES_CSF = dpwRESXp_CSF;
    SPEC.Y_pwRES_CSF = dpwRESYp_CSF;    

    PW.dt     = TempTreMethod;
    PW.pwmeth = pwdmethod;
    PW.fwhm   = lFWHM;
    PW.MAp    = MPparamNum;
    PW.ARp    = Mord;

    MatFileName = [Path2Dir '/ED' EDtype '_' num2str(BCl) '_' pwdmethod '_AR' num2str(Mord) '_MA' num2str(MPparamNum) '_FWHM' num2str(lFWHM) '_' TempTreMethod num2str(NumTmpTrend) '_ICACLEAN' num2str(icaclean) '_GSR' num2str(gsrflag) '.mat'];
    save(MatFileName,'SPEC','PW')
end

disp('xxDONExx')
