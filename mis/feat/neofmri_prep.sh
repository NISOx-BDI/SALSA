#bin/bash

set -e

SubID=$1
SesID=$2

#example: sh neofmri_prep.sh CC00698XX23 220400

#StrgDir=~/Home/GitClone/FILM2/Externals/neofmri/
StrgDir=/well/nichols/users/kfh142/data/baby/neofmri_2nd_release_rerun2

maskingflag=0
segflag=0
smoothflag=1
##################################################### Work
DataDir=${StrgDir}/sub-${SubID}/ses-${SesID}

# Segmentation file
AllSegFile=${DataDir}/fix/func_dseg

# 4D images
MCDCImgPath=${DataDir}/mcdc/func_mcdc_masked
FIXImgPath=${DataDir}/fix/func_clean_masked

cd ${DataDir}
pwd

if [ $maskingflag == 1 ]; then
	echo "neofmri_prep:: Make a brain mask from the seg file"
	# Make a brain mask from the seg files and mask the functional data
	${FSLDIR}/bin/fslmaths ${AllSegFile} -bin ${AllSegFile}_mask

	echo "neofmri_prep:: masking the MCDC and FIX images "
	${FSLDIR}/bin/fslmaths ${MCDCImgPath} -mas ${AllSegFile}_mask ${MCDCImgPath}_brain
	${FSLDIR}/bin/fslmaths ${FIXImgPath}  -mas ${AllSegFile}_mask ${FIXImgPath}_brain
fi


if [ $segflag == 1 ]; then
	# Get the WM mask isolated
	#echo "neofmri_prep:: Get the WM isolated"
	#${FSLDIR}/bin/fslmaths ${AllSegFile} -thr 3 -uthr 3 ${AllSegFile}_wm

	# Get three-tissue mask
	for i in $(seq 1 9); do
		${FSLDIR}/bin/fslmaths ${AllSegFile} -thr $i -uthr $i -bin ${DataDir}/fix/tmp_seg$i
	done

	# FAST order of tissues: CSF: 1, GM: 2, WM: 3
	# Here is the WM
	mv ${DataDir}/fix/tmp_seg3.nii.gz ${AllSegFile}_wm.nii.gz

	${FSLDIR}/bin/fslmaths ${AllSegFile}_wm -mul 3 ${AllSegFile}_wm
	# Form the graymatter
	${FSLDIR}/bin/fslmaths ${DataDir}/fix/tmp_seg2 -add ${DataDir}/fix/tmp_seg5 -add ${DataDir}/fix/tmp_seg6 -add ${DataDir}/fix/tmp_seg8 -add ${DataDir}/fix/tmp_seg9 -mul 2 ${AllSegFile}_gm
	# Form the CSF
	${FSLDIR}/bin/fslmaths ${DataDir}/fix/tmp_seg1 -add ${DataDir}/fix/tmp_seg4 -add ${DataDir}/fix/tmp_seg7 -mul 1 ${AllSegFile}_csf

	# Put all of them together
	${FSLDIR}/bin/fslmaths ${AllSegFile}_wm -add ${AllSegFile}_gm -add ${AllSegFile}_csf ${AllSegFile}_seg

	${FSLDIR}/bin/fslmaths ${AllSegFile}_wm  -bin ${AllSegFile}_wm
	${FSLDIR}/bin/fslmaths ${AllSegFile}_gm  -bin ${AllSegFile}_gm
	${FSLDIR}/bin/fslmaths ${AllSegFile}_csf -bin ${AllSegFile}_csf

	#Clean up
	/bin/rm ${DataDir}/fix/tmp_*.nii.gz
fi

if [ $smoothflag == 1 ]; then
	FWHMl=5

        origimage=${MCDCImgPath}_brain
        finalresult=${origimage}_fwhm${FWHMl}

        maskimage=${AllSegFile}_mask

        smoothparSig=$(bc -l <<< "${FWHMl}/2.3548")

        echo "===SMOOTHING:"
        echo "-IN:   ${origimage}"
        echo "-MASK: ${maskimage}"
        echo "-OUT:  ${finalresult}"
        echo "-KERNEL(mm): ${FWHMl} , sigma: ${smoothparSig}"
        echo ""

        ${FSLDIR}/bin/fslmaths $origimage -s $smoothparSig -mas $maskimage tmp_result1_tmp
        ${FSLDIR}/bin/fslmaths $maskimage -s $smoothparSig -mas $maskimage $finalresult
        ${FSLDIR}/bin/fslmaths tmp_result1_tmp -div $finalresult $finalresult

        PREFILTIMG=${finalresult}

        ${FSLDIR}/bin/imrm tmp_result1_tmp

fi
echo "  DONE"
