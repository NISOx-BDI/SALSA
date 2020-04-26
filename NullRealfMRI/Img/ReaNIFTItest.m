clear

addpath('/Users/sorooshafyouni/Home/matlab/spm12')
addpath('/Users/sorooshafyouni/Home/GitClone/FILM2/mis')
#A = '../../ExampleData/R.mpp/sub-A00029304/ses-DS2/sub-A00029304_ses-DS2_task-rest_acq-645_bold.mpp/prefiltered_func_data_bet.nii.gz';
A='/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/R_test/sub-A00027167_ses-DS2_task-rest_acq-645_bold_mpp/test.nii.gz';
system(['gunzip ' A])
%V=spm_vol(A);
%Y=spm_read_vols(V);
%YY = mean(Y,4);
[Y,Stat] = CleanNIFTI_spm('/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/R_test/sub-A00027167_ses-DS2_task-rest_acq-645_bold_mpp/test.nii');
YY = mean(Y,2);

#system(['gzip ' A]);

#V1 = Stat.spmV(1);
#V1.Removables = Stat.Removables;
#V1.fname='outputfile5.nii';

#[AA,BB] = CleanNIFTI_spm(YY,'ImgInfo',V1);
