#bin/bash
set -e

# SA,Ox,2020

######################################
#COHORT=ROCKLAND
#T=240
#TRs=0.645 #2.5 # in second
#NSUB=25
#TaskName=CHECKERBOARD
#StimulName=""
######################################
COHORT=tHCP
T=284
TRs=0.72 #2.5 # in second
NSUB=25
TaskName=MOTOR
StimulName=lh
######################################
#COHORT=tHCP
#T=253
#TRs=0.72 #2.5 # in second
#NSUB=25
#TaskName=GAMBLING
#StimulName=""
######################################

FWHMl=0
Inf2Method=FLAME1 #TFCE FLAME1

# -----------------------------------------------------------
DataDir="/well/nichols/users/scf915/${COHORT}"

TR=$(echo $TRs*1000 | bc | awk -F'.' {'print $1'})
FileNameParent=${COHORT}_${TR}_${T}_${TaskName}
ACFARO=$(echo "2*sqrt($T)" | bc)

LB=$(cat ${DataDir}/${FileNameParent}_FWHM${FWHMl}_tavg.txt)
UB=15
echo "LB: ${LB}, UB: ${UB}"

# Display ---------------------------------------------------
if [ $Inf2Method == FLAME1 ]; then
	${FSLDIR}/bin/fsleyes \
	${FSLDIR}/data/standard/MNI152_T1_2mm_brain.nii.gz \
	${DataDir}/R.PW/${FileNameParent}_gACFadjxACFadj2tT1S0P5_AR--2_MA-0_FWHM${FWHMl}_hpfk_gsr0_aroma0/LEVEL2/${Inf2Method}/tstat1.nii.gz -dr ${LB} ${UB} -cm red \
	${DataDir}/R.PW/${FileNameParent}_gFAST_AR-1_MA-0_FWHM${FWHMl}_dct_gsr0_aroma0/LEVEL2/${Inf2Method}/tstat1.nii.gz -dr ${LB} ${UB} -cm blue \
	${DataDir}/R.PW/${FileNameParent}_AR-W_AR-1_MA-0_FWHM${FWHMl}_poly_gsr0_aroma0/LEVEL2/${Inf2Method}/tstat1.nii.gz -dr ${LB} ${UB} -cm green \
	${DataDir}/R.PW/${FileNameParent}_ACF_AR-${ACFARO}_MA-0_FWHM${FWHMl}_hpf_gsr0_aroma0/LEVEL2/${Inf2Method}/tstat1.nii.gz -dr ${LB} ${UB} -cm yellow & 

elif [ $Inf2Method == TFCE ]; then
	${FSLDIR}/bin/fsleyes ${FSLDIR}/data/standard/MNI152_T1_2mm_brain.nii.gz \
	${DataDir}/R.PW/${FileNameParent}_gACFadjxACFadj2tT1S0P5_AR--2_MA-0_FWHM${FWHMl}_hpfk_gsr0_aroma0/LEVEL2/${Inf2Method}/*_Wcbhat_tstat1.nii.gz -dr ${LB} ${UB} -cm red \
	${DataDir}/R.PW/${FileNameParent}_gFAST_AR-1_MA-0_FWHM${FWHMl}_dct_gsr0_aroma0/LEVEL2/${Inf2Method}/*_Wcbhat_tstat1.nii.gz -dr ${LB} ${UB} -cm blue \
	${DataDir}/R.PW/${FileNameParent}_AR-W_AR-1_MA-0_FWHM${FWHMl}_poly_gsr0_aroma0/LEVEL2/${Inf2Method}/*_Wcbhat_tstat1.nii.gz -dr ${LB} ${UB} -cm green \
	${DataDir}/R.PW/${FileNameParent}_ACF_AR-${ACFARO}_MA-0_FWHM${FWHMl}_hpf_gsr0_aroma0/LEVEL2/${Inf2Method}/*_Wcbhat_tstat1.nii.gz -dr ${LB} ${UB} -cm yellow & 
fi
