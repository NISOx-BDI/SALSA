#bin/bash

set -e

SubID=$1
SesID=$2

#example: sh neofmri_prep.sh CC00698XX23 220400

StrgDir=~/Home/GitClone/FILM2/Externals/neofmri/

##################################################### Work
DataDir=${StrgDir}/sub-${SubID}/ses-${SesID}

# Segmentation file
AllSegFile=${DataDir}/fix/func_dseg

# 4D images
MCDCImgPath=${DataDir}/mcdc/func_mcdc_masked
FIXImgPath=${DataDir}/fix/func_clean_masked

echo "neofmri_prep:: Make a brain mask from the seg file"
# Make a brain mask from the seg files and mask the functional data
${FSLDIR}/bin/fslmaths ${AllSegFile} -bin ${AllSegFile}_mask

echo "neofmri_prep:: masking the MCDC and FIX images "
${FSLDIR}/bin/fslmaths ${MCDCImgPath} -mas ${AllSegFile}_mask ${MCDCImgPath}_brain
${FSLDIR}/bin/fslmaths ${FIXImgPath}  -mas ${AllSegFile}_mask ${FIXImgPath}_brain

# Get the WM mask isolated
echo "neofmri_prep:: Get the WM isolated"
${FSLDIR}/bin/fslmaths ${AllSegFile} -thr 3 -uthr 3 ${AllSegFile}_wm
