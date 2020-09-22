#bin/bash

######################################
COHORT=ROCKLAND
T=240
TRs=0.645 #2.5 # in second
NSUB=25
TaskName=CHECKERBOARD
StimulName=""
######################################
#COHORT=tHCP
#T=284
#TRs=0.72 #2.5 # in second
#NSUB=25
#TaskName=MOTOR
#StimulName=lh
######################################
#COHORT=tHCP
#T=253
#TRs=0.72 #2.5 # in second
#NSUB=25
#TaskName=GAMBLING
#StimulName=""
######################################

FWHMl=5

DataDir="/well/nichols/users/scf915/${COHORT}"

TR=$(echo $TRs*1000 | bc | awk -F'.' {'print $1'})
FileNameParent=${COHORT}_${TR}_${T}_${TaskName}
ACFARO=$(echo "2*sqrt($T)" | bc)

cat ${DataDir}/R.PW/${FileNameParent}_*_AR-*_MA-0_FWHM${FWHMl}_*_gsr0_aroma0/LEVEL2/TFCE/*.log | grep "Wcbhat_vox_tstat1 is:" | awk '{ print $6 }' \
>> ${DataDir}/${FileNameParent}_FWHM${FWHMl}_gsr0_aroma0_TFCE_tstat_criticalval.txt

cat ${DataDir}/${FileNameParent}_FWHM${FWHMl}_gsr0_aroma0_TFCE_tstat_criticalval.txt | awk '{ total += $1; count++ } END { print total/count }' \
> ${DataDir}/${FileNameParent}_FWHM${FWHMl}_tavg.txt

