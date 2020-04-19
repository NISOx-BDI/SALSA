#!bin/bash
#
# Calculates the bias of the standard error of the Bhat
# Bias=(SE-SIM/SIM)
#
#



set -e

T=200
TR=645
COHORT=ROCKLAND

fwhml=0

pwmethodlist=(AR-W AR-YW ACF ARMAHR)

for pwmethod in ${pwmethodlist[@]}
do

	MAO=0
	[ $pwmethod == ARMAHR ]&& MAO=1

	AROlist=(1 2 5 10 20)
	[ $pwmethod == ACF ]&& AROlist=(5 10 15 30 60)

	for ARO in ${AROlist[@]}
	do
		echo "${COHORT}_${T}_${TR} ${pwmethod} AR-${ARO} MA-${MAO} FWHM${fwhml}"

		SIMDIR=/well/nichols/users/scf915/${COHORT}/R.SIM/SIM_${COHORT}_${T}_${TR}_${pwmethod}_AR-${ARO}_MA-${MAO}_FWHM${fwhml}

#		NIMG=$(ls $Bhat_WDC | wc -l)
#		echo "Merging ${NIMG} images and summarising mean(s.e. Bhat) and std(Bhat)"

		# Get standard error of simulated Betas
		Bhat_WDC=${SIMDIR}/NullImg_*/Sim*_EDboxcar_20_${pwmethod}_AR${ARO}_MA${MAO}_Bhat_PW.nii.gz
		fslmerge -t ${SIMDIR}/Bhat_PW.nii.gz ${Bhat_WDC}
		fslmaths ${SIMDIR}/Bhat_PW.nii.gz -Tstd ${SIMDIR}/STD_Bhat_PW.nii.gz

		# Get mean of textbook betas
		seBhat_WDC=${SIMDIR}/NullImg_*/Sim*_EDboxcar_20_${pwmethod}_AR${ARO}_MA${MAO}_SE_PW.nii.gz
		fslmerge -t ${SIMDIR}/seBhat_PW.nii.gz ${seBhat_WDC}
		fslmaths ${SIMDIR}/seBhat_PW.nii.gz -Tmean ${SIMDIR}/mse_Bhat_PW.nii.gz

		# ( MeanSEBeta - STDBeta ) / STDBeta
		BiasFileName=${SIMDIR}/Bias_BhatSE_${T}_${TR}_${pwmethod}_AR${ARO}_MA${MAO}_FWHM${fwhml}.nii.gz
		fslmaths ${SIMDIR}/mse_Bhat_PW.nii.gz -sub ${SIMDIR}/STD_Bhat_PW.nii.gz -div ${SIMDIR}/STD_Bhat_PW.nii.gz ${BiasFileName}

		echo "Bias: ${BiasFileName}"
	done
done

echo "-- DONE --"
