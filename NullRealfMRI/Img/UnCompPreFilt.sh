
COHORT=ROCKLAND
SesID=DS2

SmthMM=5

RawData=/well/nichols/users/scf915/${COHORT}/R_mpp
SubIDDir=/well/nichols/users/scf915/${COHORT}


module load fsl
source ${HOME}/.fslconf/fsl.sh

cat ${SubIDDir}/sub.txt | head -n 200 | while read SubID
do
	echo ${SubID}
	SubID=$(echo ${SubID} | awk -F"-" '{print $2}')
	echo ${SubID}

	DataDir=${RawData}/sub-${SubID}/ses-${SesID}/sub-${SubID}_ses-DS2_task-rest_acq-645_bold_mpp/
	CompData=${DataDir}/prefiltered_func_data_bet

	#echo ${CompData}

	if [ ! -f ${CompData}.nii ]; then
		echo "Uncompressing..."
		gunzip -c ${CompData}.nii.gz #> ${CompData}.nii
	elif [ -f ${CompData}.nii ] && [ -f ${CompData}.nii.gz ]; then
		rm ${CompData}.nii.gz
	fi

#	if [ ! -f ${CompData}_FWHM${SmthMM}.nii ]; then

		${FSLDIR}/bin/fslmaths ${CompData}.nii -kernel gauss 2.1233226 -fmean ${CompData}_FWHM${SmthMM}.nii
#		echo "SMOOTHED AND DECOMPRESSED: ${CompData}_FWHM${SmthMM}.nii"
#	elif [ ${CompData}_FWHM${SmthMM}.nii.gz ]; then
#		rm -f ${CompData}_FWHM${SmthMM}.nii.gz
#	fi

done


