
#-Wiktor Olszowy

#-simulate 4D fMRI scans similar to 'FCP_Beijing', which is resting-state data

library(neuRosim)
library(oro.nifti)
# library(parallel)
# library(AnalyzeFMRI)

setwd("/Users/sorooshafyouni/Home/PREW/SimulatefMRI")

FCP_Beijing_RS = readNIfTI("/Users/sorooshafyouni/Home/PREW/SimulatefMRI/R/sub-A00029304/ses-DS2/func/sub-A00029304_ses-DS2_task-rest_acq-1400_bold.nii.gz")
baseline       = apply(FCP_Beijing_RS, 1:3, mean)
abbr           = "SIM"

#res = mclapply(1:100, function(sub_id) {
  sub_id=1
  TRVal=2
  Nt=225
   #-the following ordering, taken from manual, is wrong!
   #-weights refer to 6 different noises: "white", "temporal", "spatial", "low-frequency", "physiological", "task-related"
   #-from code:
   #-n <- (w[1] * n.white + w[2] * n.temp + w[3] * n.low + w[4] * n.phys + w[5] * n.task + w[6] * n.spat)/sqrt(sum(w^2))
   A = simVOLfmri(dim=dim(baseline), nscan=Nt, TR=TRVal, noise="mixture", rho.temp=0.48,
      w=c(0.25, 0.50, 0, 0, 0, 0.25), SNR=10, base=baseline)
   #-so that later ordering is easy
   sub_id_4digit = paste0('TR',TRVal*1000,'T',Nt,'_', sub_id)
   #sub_id_4digit = substr(sub_id_4digit, nchar(sub_id_4digit)-3, nchar(sub_id_4digit))
   writeNIfTI(nifti(A, datatype=16, pixdim=c(0,2,2,2,TRVal,1,1,1)), paste0("sub-", abbr, sub_id_4digit, "_rest_bold"))
   system(paste0("gunzip sub-", abbr, sub_id_4digit, "_rest_bold.nii.gz"))
   #system(paste0("cp Beijing_sub98617_T1.nii sub-",       abbr, sub_id_4digit, "_T1w.nii"))
   #system(paste0("cp Beijing_sub98617_T1_brain.nii sub-", abbr, sub_id_4digit, "_T1w_brain.nii"))

#}, mc.cores=24)

#A = f.read.nifti.volume("Beijing_sub98617.nii")
#A = f.read.nifti.volume("sub-SIM0001_rest_bold.nii")
#dim(A)
#all_AR1 = c(0)
#for (i in 20:25) {
#   for (j in 20:25) {
#      for (k in 10:14) {
#         ts = A[i,j,k,]
#         #cat(arima(ts, order=c(1,0,0))$coef[1], "\n")
#         all_AR1 = c(all_AR1, arima(ts, order=c(1,0,0))$coef[1])
#      }
#   }
#}
#mean(all_AR1)
