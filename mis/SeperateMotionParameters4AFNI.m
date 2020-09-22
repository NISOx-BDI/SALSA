clear

COHORT = 'ROCKLAND'; 
TR        = 0.645;
MParamNum = 24;
COHORTDIR = ['/well/nichols/users/scf915/' COHORT];

disp('=====Set Path=============================')
PATH2AUX='/Users/sorooshafyouni/Home/GitClone/FILM2';
addpath([PATH2AUX '/mis'])

disp('=====Run Stuff=============================')

load([COHORTDIR '/' COHORT '_sesid.mat'])
load([COHORTDIR '/' COHORT '_subid.mat'])
SubList = participants;
SesList = ses;

for subcnt = 1:60
    
    SubID = SubList{subcnt};
    SesID = SesList{subcnt};
    
    %Raw Images (MMP feat output)
    Path2ImgRaw = [COHORTDIR '/R_mpp/sub-' SubID '/ses-' SesID];
    if strcmpi(COHORT,'ROCKLAND')
        Path2ImgDir = [Path2ImgRaw '/sub-' SubID '_ses-' SesID '_task-rest_acq-' num2str(TR*1000) '_bold_mpp'];
    elseif any(strcmpi(COHORT,{'Beijing','Cambridge'}))
        Path2ImgDir = [Path2ImgRaw '/rest_mpp'];
    end

    Path2MC   = [Path2ImgDir '/prefiltered_func_data_mcf.par'];

    disp(['Image: ' Path2Img])
    disp(['Motion params: ' Path2MC])

    % Motion parameters ----------------------------------------

    MCp = load(Path2MC);
    MCp = GenMotionParam(MCp,MParamNum); 
    
    path2saveMVs = [Path2ImgRaw '/mp'];
    if ~exist(path2saveMVs, 'dir')
        mkdir(path2saveMVs)
    end
    
    for mvcnt=1:MParamNum
        filename=[path2saveMVs '/' COHORT '_SubID-' SubID '_ses-' SesID '_T' num2str(T) '_TR' num2str(TR*1000) '_' num2str(mvcnt) '.txt'];
        fileID = fopen(filename,'w');
        fprintf(fileID,'%f\n',MCp(:,mvcnt));
        fclose(fileID);        
    end

end