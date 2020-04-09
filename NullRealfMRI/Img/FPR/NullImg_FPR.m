clear


get_tcrit = @(alph,df) tinv(1-alph,df);

T       = 900;
rankX   = 35;
alphlev = 0.05;

df      = T-rankX;
t_crit  = get_tcrit(alphlev,df);


Y = CleanNIFTI_spm('/Users/sorooshafyouni/Home/GitClone/FILM2/ExampleData/R.mpp/RNullfMRI_A00029304_DS2/AR-W_AR5_tVALUE_Naive.nii');

num_sigvox  = sum(abs(Y)>t_crit);
FPR         = num_sigvox./numel(Y); 