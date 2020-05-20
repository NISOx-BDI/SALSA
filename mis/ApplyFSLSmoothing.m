function [sY,smStat] = ApplyFSLSmoothing(Y,FWHMl,ImgStat,path2mask)
    
    % where is the fslsmooth.sh code
    path2sh = which('ApplyFSLSmoothing');
    spathstr=fileparts(path2sh);
    
    % Make FSL happy
    %system(['FSLOUTPUTTYPE=NIFTI_PAIR; export FSLOUTPUTTYPE; /bin/sh /${FSLDIR}/etc/fslconf/fsl.sh']);
    
    %  make a temp nifti file out of the data
    CLK     = fix(clock);
    tmpdir  = [tempdir 'octspm12/tmp_' num2str(randi(5000)) '_' num2str(CLK(end))]; 
    mkdir(tmpdir)
    disp(['Unzip into: ' tmpdir ])
    randtempfilename2smooth=[tmpdir '/ToBesmooth_tmp_' num2str(randi(50)) num2str(CLK(end)+randi(10)) '.nii'];
    randtempfilenamesmooth=[tmpdir '/smooth_tmp_' num2str(randi(50)) num2str(CLK(end)+randi(10)) '.nii'];

    % build back a nifti file
    CleanNIFTI_spm(Y,'ImgInfo',ImgStat.spmV,'destdir',randtempfilename2smooth,'Removables',ImgStat.Removables,'verbose',0);
    
    fslsmoothcmd = ['sh ' spathstr '/fslsmooth.sh ' randtempfilename2smooth ' ' path2mask ' ' num2str(FWHMl) ' ' randtempfilenamesmooth];
    
    % Do the smoothing
    disp(fslsmoothcmd)
    status = call_fsl(fslsmoothcmd);
    
    if ~status
        system(['gunzip ' randtempfilenamesmooth]);
        [sY,smStat] = CleanNIFTI_spm(randtempfilenamesmooth,'mask',path2mask,'verbose',0);
        system(['rm -r ' tmpdir]);
    elseif status
        error('ApplyFSLSmoothing:: fsl_call failed.')
    end
    
end


% CLK	 = fix(clock);
% tmpdir  = [tempdir 'octspm12/tmp_' num2str(randi(5000)) '_' num2str(CLK(end))]; % make a temp directory 
% mkdir(tmpdir)
% disp(['Unzip into: ' tmpdir ])
% randtempfilename=[tmpdir '/prefilt_tmp_' SubID '_' num2str(randi(50)) num2str(CLK(end)+randi(10)) '.nii'];
% system(['gunzip -c ' Path2Img ' > ' randtempfilename]); %gunzip function in Octave deletes the source file in my version!
% 
% [Y,InputImgStat]=CleanNIFTI_spm(randtempfilename,'demean','fwhm',lFWHM);
% 
% disp(['Remove the temp directory: ' tmpdir])
% %rmdir(tmpdir,'s')
% system(['rm -rf ' tmpdir])