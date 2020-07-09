#bin/bash

set -e

SubID=$1
SesID=$2

#example: sh neofmri_prep.sh CC00698XX23 220400

StrgDir=~/Home/GitClone/FILM2/Externals/neofmri/

maskingflag=0
##################################################### Work
DataDir=${StrgDir}/sub-${SubID}/ses-${SesID}

# Segmentation file
AllSegFile=${DataDir}/fix/func_dseg

# 4D images
MCDCImgPath=${DataDir}/mcdc/func_mcdc_masked
FIXImgPath=${DataDir}/fix/func_clean_masked

if [ $maskingflag == 1 ]; then
	echo "neofmri_prep:: Make a brain mask from the seg file"
	# Make a brain mask from the seg files and mask the functional data
	${FSLDIR}/bin/fslmaths ${AllSegFile} -bin ${AllSegFile}_mask

	echo "neofmri_prep:: masking the MCDC and FIX images "
	${FSLDIR}/bin/fslmaths ${MCDCImgPath} -mas ${AllSegFile}_mask ${MCDCImgPath}_brain
	${FSLDIR}/bin/fslmaths ${FIXImgPath}  -mas ${AllSegFile}_mask ${FIXImgPath}_brain
fi

# Get the WM mask isolated
#echo "neofmri_prep:: Get the WM isolated"
#${FSLDIR}/bin/fslmaths ${AllSegFile} -thr 3 -uthr 3 ${AllSegFile}_wm

# Get three-tissue mask
for i in $(seq 1 9); do
	fslmaths ${AllSegFile} -thr $i -uthr $i -bin ${DataDir}/fix/tmp_seg$i
done

# FAST order of tissues: CSF: 1, GM: 2, WM: 3
# Here is the WM
mv ${DataDir}/fix/tmp_seg3.nii.gz ${AllSegFile}_wm.nii.gz

fslmaths ${AllSegFile}_wm -mul 3 ${AllSegFile}_wm
# Form the graymatter
fslmaths ${DataDir}/fix/tmp_seg2 -add ${DataDir}/fix/tmp_seg5 -add ${DataDir}/fix/tmp_seg6 -add ${DataDir}/fix/tmp_seg8 -add ${DataDir}/fix/tmp_seg9 -mul 2 ${AllSegFile}_gm
# Form the CSF
fslmaths ${DataDir}/fix/tmp_seg1 -add ${DataDir}/fix/tmp_seg4 -add ${DataDir}/fix/tmp_seg7 -mul 1 ${AllSegFile}_csf

# Put all of them together
fslmaths ${AllSegFile}_wm -add ${AllSegFile}_gm -add ${AllSegFile}_csf ${AllSegFile}_seg

