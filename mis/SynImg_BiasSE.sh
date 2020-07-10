#!bin/bash
#
# Calculates the bias of the standard error of the Bhat
# Bias=(SE-SIM/SIM)
#
#

module load fsl

set -e

T=900
TR=645
COHORT=ROCKLAND

fwhml=0

TempTreMethod="hpf"
NumTempTre=100

pwmethodlist=(ACFadj ARMAHR) #AR-W ACF)

for pwmethod in ${pwmethodlist[@]}
do

	MAO=0
	[ $pwmethod == ARMAHR ]&& MAO=1

	AROlist=(1 2 5 10 20)
	[ $pwmethod == ACF ] || [ $pwmethod == ACFadj ]&& AROlist=(5 10 15 $(echo "sqrt($T)" | bc) $(echo "2*sqrt($T)" | bc))

	for ARO in ${AROlist[@]}
	do
		echo "${COHORT}_${T}_${TR} ${pwmethod} AR-${ARO} MA-${MAO} FWHM${fwhml} ${TempTreMethod}${NumTempTre}"

		SIMDIR=/well/nichols/users/scf915/${COHORT}/R.SIM/SIM_${COHORT}_${T}_${TR}_${pwmethod}_AR-${ARO}_MA-${MAO}_FWHM${fwhml}_${TempTreMethod}

#		NIMG=$(ls $Bhat_WDC | wc -l)
#		echo "Merging ${NIMG} images and summarising mean(s.e. Bhat) and std(Bhat)"

		# Get standard error of simulated Betas
		Bhat_WDC=${SIMDIR}/NullImg_*/EDboxcar_20_${pwmethod}_AR${ARO}_MA${MAO}_FWHM0_${TempTreMethod}${NumTempTre}_Bhat_PW.nii.gz
		${FSLDIR}/bin/fslmerge -t ${SIMDIR}/Bhat_PW.nii.gz ${Bhat_WDC}
		${FSLDIR}/bin/fslmaths ${SIMDIR}/Bhat_PW.nii.gz -Tstd ${SIMDIR}/STD_Bhat_PW.nii.gz

		echo "----- Merger of Bhat_PW is done."

		# Get mean of textbook betas
		#seBhat_WDC=${SIMDIR}/NullImg_*/Sim*_EDboxcar_20_${pwmethod}_AR${ARO}_MA${MAO}_SE_PW.nii.gz
		#_FWHM0_hpf100_CPZ_PW.nii.gz
		seBhat_WDC=${SIMDIR}/NullImg_*/EDboxcar_20_${pwmethod}_AR${ARO}_MA${MAO}_FWHM0_${TempTreMethod}${NumTempTre}_SE_PW.nii.gz
		${FSLDIR}/bin/fslmerge -t ${SIMDIR}/seBhat_PW.nii.gz ${seBhat_WDC}
		${FSLDIR}/bin/fslmaths ${SIMDIR}/seBhat_PW.nii.gz -Tmean ${SIMDIR}/mse_Bhat_PW.nii.gz

		echo "----- Merger of SE_PW is done."

		# ( MeanSEBeta - STDBeta ) / STDBeta
		BiasFileName=${SIMDIR}/Bias_BhatSE_${T}_${TR}_${pwmethod}_AR${ARO}_MA${MAO}_FWHM${fwhml}_${TempTreMethod}${NumTempTre}.nii.gz
		${FSLDIR}/bin/fslmaths ${SIMDIR}/mse_Bhat_PW.nii.gz -sub ${SIMDIR}/STD_Bhat_PW.nii.gz -div ${SIMDIR}/STD_Bhat_PW.nii.gz ${BiasFileName}

		echo "Bias: ${BiasFileName}"

		tVALUEFILENAME_PW=${SIMDIR}/tVALUE_PW_${T}_${TR}_${pwmethod}_AR${ARO}_MA${MAO}_FWHM${fwhml}_${TempTreMethod}${NumTempTre}.nii.gz
		${FSLDIR}/bin/fslmerge -t ${tVALUEFILENAME_PW} ${SIMDIR}/NullImg_*/EDboxcar_20_${pwmethod}_AR${ARO}_MA${MAO}_FWHM0_${TempTreMethod}${NumTempTre}_tVALUE_PW.nii.gz

		echo "----- Merger of tvalue of PW is done."

		tVALUEFILENAME_Naive=${SIMDIR}/tVALUE_Naive_${T}_${TR}_${pwmethod}_AR${ARO}_MA${MAO}_FWHM${fwhml}_${TempTreMethod}${NumTempTre}.nii.gz
		${FSLDIR}/bin/fslmerge -t ${tVALUEFILENAME_Naive} ${SIMDIR}/NullImg_*/EDboxcar_20_${pwmethod}_AR${ARO}_MA${MAO}_FWHM0_${TempTreMethod}${NumTempTre}_tVALUE_Naive.nii.gz

		echo "----- Merger of tvalue of Naive is done."

	done
done

echo "-- DONE --"
