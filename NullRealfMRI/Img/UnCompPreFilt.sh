
COHORT=ROCKLAND
SesID=DS2

RawData=/well/nichols/users/scf915/${COHORT}/R_mpp
SubIDDir=/well/nichols/users/scf915/${COHORT}

cat ${SubIDDir}/sub.txt | head -n 200 | while read SubID
do
	echo ${SubID}
	SubID=$(echo ${SubID} | awk -F"-" '{print $2}')
	echo ${SubID}

	CompData=${RawData}/sub-${SubID}/ses-${SesID}/sub-${SubID}_ses-DS2_task-rest_acq-645_bold_mpp/prefiltered_func_data_bet.nii.gz
	UnCompData=${RawData}/sub-${SubID}/ses-${SesID}/sub-${SubID}_ses-DS2_task-rest_acq-645_bold_mpp/prefiltered_func_data_bet.nii

	echo ${CompData}

	if [ ! -f $UnCompData ]; then
		echo "Uncompressing..."
		gunzip -c ${CompData} > ${UnCompData}
	fi

done


