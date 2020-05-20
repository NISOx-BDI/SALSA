function [sY,smStat] = ApplyFSLSusan(Y,FWHMl,ImgStat,path2mask)
    
%     featlib.tcl:
%     set int_2_98 [ fsl:exec "${FSLDIR}/bin/fslstats $funcdata -p 2 -p 98" ]
%     set int2  [ lindex $int_2_98 0 ]
%     set median_intensity [ fsl:exec "${FSLDIR}/bin/fslstats $funcdata_unmasked -k mask -p 50" ]
% 
%     set smoothsigma [ expr $fmri(smooth) / 2.355 ]
%     set susan_int [ expr ( $median_intensity - $int2 ) * 0.75 ]
%     fsl:exec "${FSLDIR}/bin/fslmaths $funcdata -Tmean mean_func"
%     fsl:exec "${FSLDIR}/bin/susan $funcdata $susan_int $smoothsigma 3 1 1 mean_func $susan_int prefiltered_func_data_smooth"



    % where is the fslsmooth.sh code
    path2sh = which('ApplyFSLSusan');
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
    CleanNIFTI_spm(Y,'ImgInfo',ImgStat.spmV,'destdir',randtempfilename2smooth,'Removables',ImgStat.Removables);
    
    % run internal function fslsusan.sh
    fslsusancmd = ['sh ' spathstr '/fslsusan.sh ' randtempfilename2smooth ' ' path2mask ' ' num2str(FWHMl) ' ' randtempfilenamesmooth];
    
    % Do the smoothing
    disp(fslsusancmd)
    status = call_fsl(fslsusancmd);
    
    if ~status
        system(['gunzip ' randtempfilenamesmooth]);
        [sY,smStat] = CleanNIFTI_spm(randtempfilenamesmooth,'mask',path2mask);
        system(['rm -r ' tmpdir]);
    elseif status
        error('ApplyFSLSUSAN:: fsl_call failed.')
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