#bin/bash

# Inputs
SubID=$1
SesID=$2

roithr=2

AtlasMNI=${HOME}/Home/Atlas/Yeo/Yeo2011_17Networks_FSL_MNI152_2mm.nii.gz

Path2ImgDir=${HOME}/Home/GitClone/FILM2/Externals/ROCKLAND
ImgDir=${Path2ImgDir}/sub-${SubID}/ses-${SesID}/sub-${SubID}_ses-${SesID}_task-rest_acq-645_bold_mpp
RegDir=${ImgDir}/reg
#---------------------------------------------------------------------------
warpvol=${ImgDir}/reg/example_func2standard_warp.nii.gz
invwarpvol=${ImgDir}/reg/standard2example_func_warp.nii.gz
invflirtmat=${ImgDir}/reg/standard2example_func.mat
refvol=${ImgDir}/reg/example_func.nii.gz
invol=${ImgDir}/atlas/YeoMNI.nii.gz
outvol=${ImgDir}/atlas/Yeo_func

#copy atlas ---------------------------------------------------------------------------
mkdir -p ${ImgDir}/atlas/
cp $AtlasMNI $invol

# Inverse Warp ---------------------------------------------------------------------------
${FSLDIR}/bin/invwarp -w $warpvol -o $invwarpvol -r $refvol

# Inverse mat, Or does it already exists?
#${FSLDIR}/bin/convert_xfm -omat B2A.mat -inverse A2B.mat

# bring the atlas into the EPI space --------------------------------------------------------------------------
${FSLDIR}/bin/applywarp -i $invol -o $outvol -r $refvol -w $invwarpvol --interp=nn #--premat=$invflirtmat --interp=nn

# Get ROI out ---------------------------------------------------------------------------
#uthr=$(( roithr+1 ))
#lthr=$(( roithr-1 ))
#echo $lthr $uthr
${FSLDIR}/bin/fslmaths ${outvol} -thr $roithr -uthr $roithr -bin ${outvol}_$roithr
