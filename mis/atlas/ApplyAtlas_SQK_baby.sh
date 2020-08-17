#bin/bash

set -e

############################################################################
#FUNCTIONS #################################################################

function GetAtlas () {

# Inputs
Cohort=$1
SubID=$2
SesID=$3

ImgDir=/well/nichols/users/kfh142/data/baby/neofmri_2nd_release_rerun2/sub-${SubID}/ses-${SesID}/
#---------------------------------------------------------------------------
invwarpvol=${ImgDir}/reg/func-mcdc2standard_invwarp.nii.gz
refvol=${ImgDir}/ica/mean.nii.gz

#AffixList=(squaremask squaremaskkernel)
#InterList=(nn trilinear)

AffixList=(squaremaskkernel_neo)
InterList=(trilinear)

aff_cnt=0
for Affix in ${AffixList[@]}; do

	InterpTmp=${InterList[$aff_cnt]}

	echo "On: ${Affix}"
	echo "Interp: ${InterpTmp}"

	#AtlasMNI=${HOME}/bin/FILM2/mis/atlas/${Affix}.nii.gz
	#invol=${ImgDir}/atlas/${Affix}MNI.nii.gz

	invol=${HOME}/bin/FILM2/mis/atlas/${Affix}.nii.gz
	outvol=${ImgDir}/atlas/${Affix}_func

	#copy atlas ---------------------------------------------------------------------------
	mkdir -p ${ImgDir}/atlas/

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

#TR=$(echo $TRs*1000 | bc | awk -F'.' {'print $1'})

for (( s=1; s<=60; s++ ))
do
        SubID=$(cat /well/nichols/users/scf915/${COHORT}/${COHORT}_subid.txt | awk {'print $1'} | sed "${s}q;d")
	SesID=$(cat /well/nichols/users/scf915/${COHORT}/${COHORT}_sesid.txt | awk {'print $1'} | sed "${s}q;d")
	echo $COHORT $TRs $SubID $SesID

	GetAtlas $COHORT $SubID $SesID

done



