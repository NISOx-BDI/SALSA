clear

addpath('~/bin/FILM2/mis')
addpath('/well/nichols/users/scf915/externals/spm12')

#COHORTDIR='/well/nichols/users/scf915/ROCKLAND';
#Path2Img = [COHORTDIR '/R_mpp/sub-A00028994/ses-DS2/sub-A00028994_ses-DS2_task-rest_acq-1400_bold_mpp/prefiltered_func_data_bet.nii.gz'];

Path2Img='/well/nichols/users/scf915/rawHCP/943862_3T_rfMRI_REST1_LR.nii.gz';


CLK	 = fix(clock);
tmpdir   = [tempdir '/octspm12/tmp_' num2str(randi(5000)) '_' num2str(CLK(end))]; % make a temp directory 
mkdir(tmpdir)
randtempfilename = [tmpdir '/prefilt_tmp_' num2str(randi(50)) num2str(CLK(end)+randi(10)) '.nii'];
disp(['Unzip into: ' randtempfilename ])
system(['gunzip -c ' Path2Img ' > ' randtempfilename]); %gunzip function in Octave deletes the source file in my version!
[Y,InputImgStat] = CleanNIFTI_spm(randtempfilename,'demean');
disp(['Remove the temp directory: ' tmpdir])
system(['rm -rf ' tmpdir])
%rmdir(tmpdir,'s')


%disp(A)
%[Y,Stat] = CleanNIFTI_spm(A);
%YY = mean(Y,2);

disp('DONE')

#V1 = Stat.spmV(1);
#V1.Removables = Stat.Removables;
#V1.fname='outputfile5.nii';
#[AA,BB] = CleanNIFTI_spm(YY,'ImgInfo',V1);
