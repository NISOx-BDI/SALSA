#bin/bash



############################################################################
#FUNCTIONS #################################################################

function GetAtlas () {

# Inputs
Cohort=$1
TRs=$2
SubID=$3
SesID=$4

roithr=2


TR=$(echo $TRs*1000 | bc | awk -F'.' {'print $1'})

#echo $TR

AtlasMNI=${HOME}/bin/FILM2/mis/atlas/Yeo2011_17Networks_FSL_MNI152_2mm.nii.gz

Path2ImgDir=/well/nichols/users/scf915/${Cohort}/R_mpp

if [ $Cohort == ROCKLAND ]; then
	ImgDir=${Path2ImgDir}/sub-${SubID}/ses-${SesID}/sub-${SubID}_ses-${SesID}_task-rest_acq-${TR}_bold_mpp
elif [ $Cohort == Beijing ]; then
	ImgDir=${Path2ImgDir}/sub-${SubID}/ses-${SesID}/rest_mpp
elif [ $Cohort == HCP ]; then
	ImgDir=${Path2ImgDir}/sub-${SubID}/ses-${SesID}/${SubID}_3T_rfMRI_${SesID}_mpp
fi


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

}


#MAIN ######################################################################
############################################################################

COHORT=$1
TRs=$2
SesID=$3

#TR=$(echo $TRs*1000 | bc | awk -F'.' {'print $1'})

for (( s=1; s<=51; s++ ))
do
        SubID=$(cat /well/nichols/users/scf915/${COHORT}/${COHORT}_subid.txt | awk {'print $1'} | sed "${s}q;d")

	echo $COHORT $TRs $SubID $SesID

	GetAtlas $COHORT $TRs $SubID $SesID

done



