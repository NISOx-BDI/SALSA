#bin/bash


############################################################################
#FUNCTIONS #################################################################

function GetAtlas () {

# Inputs
Cohort=$1
TRs=$2
SubID=$3
SesID=$4


if [ $Cohort == ROCKLAND ]; then
	Path2ImgDir=/well/nichols/users/scf915/${Cohort}/R_mpp
	#TRs=0.645
	TR=$(echo $TRs*1000 | bc | awk -F'.' {'print $1'})
	ImgDir=${Path2ImgDir}/sub-${SubID}/ses-${SesID}/sub-${SubID}_ses-${SesID}_task-rest_acq-${TR}_bold_mpp
elif [ $Cohort == Beijing ]; then
	Path2ImgDir=/well/nichols/users/scf915/${Cohort}/R_mpp
	#TRs=
	TR=$(echo $TRs*1000 | bc | awk -F'.' {'print $1'})
	ImgDir=${Path2ImgDir}/sub-${SubID}/ses-${SesID}/rest_mpp
elif [ $Cohort == HCP ]; then
	Path2ImgDir=/well/nichols/users/scf915/${Cohort}/R_mpp
	#TRs=
	TR=$(echo $TRs*1000 | bc | awk -F'.' {'print $1'})
	ImgDir=${Path2ImgDir}/sub-${SubID}/ses-${SesID}/${SubID}_3T_rfMRI_${SesID}_mpp
fi


RegDir=${ImgDir}/reg
#---------------------------------------------------------------------------
warpvol=${ImgDir}/reg/example_func2standard_warp.nii.gz
invwarpvol=${ImgDir}/reg/standard2example_func_warp.nii.gz
invflirtmat=${ImgDir}/reg/standard2example_func.mat
refvol=${ImgDir}/reg/example_func.nii.gz

#AffixList=(squaremask squaremaskkernel)
#InterList=(nn trilinear)

AffixList=(squaremaskkernel)
InterList=(trilinear)

aff_cnt=0
for Affix in ${AffixList[@]}; do

	InterpTmp=${InterList[$aff_cnt]}

	echo "On: ${Affix}"
	echo "Interp: ${InterpTmp}"

	AtlasMNI=${HOME}/bin/FILM2/mis/atlas/${Affix}.nii.gz
	invol=${ImgDir}/atlas/${Affix}MNI.nii.gz
	outvol=${ImgDir}/atlas/${Affix}_func

	#copy atlas ---------------------------------------------------------------------------
	mkdir -p ${ImgDir}/atlas/
	cp $AtlasMNI $invol

	# Inverse Warp ---------------------------------------------------------------------------
	${FSLDIR}/bin/invwarp -w $warpvol -o $invwarpvol -r $refvol

	# Inverse mat, Or does it already exists?
	#${FSLDIR}/bin/convert_xfm -omat B2A.mat -inverse A2B.mat

	# bring the atlas into the EPI space --------------------------------------------------------------------------
	${FSLDIR}/bin/applywarp -i $invol -o $outvol -r $refvol -w $invwarpvol --interp=$InterpTmp

	aff_cnt=$((aff_cnt+1 ))

	# Get ROI out ---------------------------------------------------------------------------
	#uthr=$(( roithr+1 ))
	#lthr=$(( roithr-1 ))
	#echo $lthr $uthr
	${FSLDIR}/bin/fslmaths ${outvol} -bin ${outvol}_mask
done

}


#MAIN ######################################################################
############################################################################

COHORT=$1
TRs=$2

#TR=$(echo $TRs*1000 | bc | awk -F'.' {'print $1'})

for (( s=1; s<=51; s++ ))
do
        SubID=$(cat /well/nichols/users/scf915/${COHORT}/${COHORT}_subid.txt | awk {'print $1'} | sed "${s}q;d")
	SesID=$(cat /well/nichols/users/scf915/${COHORT}/${COHORT}_sesid.txt | awk {'print $1'} | sed "${s}q;d")
	echo $COHORT $TRs $SubID $SesID

	GetAtlas $COHORT $TRs $SubID $SesID

done



