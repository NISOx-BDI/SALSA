clear

addpath('~/bin/FILM2/mis')
addpath('/well/nichols/users/scf915/externals/spm12')

COHORTDIR='/well/nichols/users/scf915/ROCKLAND';
A = [COHORTDIR '/R_mpp/sub-A00029304/ses-DS2/sub-A00029304_ses-DS2_task-rest_acq-645_bold_mpp/prefiltered_func_data_bet.nii'];

disp(A)

[Y,Stat] = CleanNIFTI_spm(A);
YY = mean(Y,2);

#V1 = Stat.spmV(1);
#V1.Removables = Stat.Removables;
#V1.fname='outputfile5.nii';
#[AA,BB] = CleanNIFTI_spm(YY,'ImgInfo',V1);
