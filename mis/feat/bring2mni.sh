


function ToStandard () {

invol=$1

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
refvol=${ImgDir}/reg/example_func.nii.gz


invol=${ImgDir}/atlas/YeoMNI.nii.gz
outvol=${ImgDir}/atlas/Yeo_func

#copy atlas ---------------------------------------------------------------------------
mkdir -p ${ImgDir}/mni/

# bring the atlas into the EPI space --------------------------------------------------------------------------
${FSLDIR}/bin/applywarp -i $invol -o $outvol -r $refvol -w $warpvol #--interp=nn #--premat=$invflirtmat --interp=nn

}



#MAIN ######################################################################
############################################################################

COHORT=$1
TRs=$2
SesID=$3

#TR=$(echo $TRs*1000 | bc | awk -F'.' {'print $1'})

for (( s=1; s<=51; s++ ))
do
        #SubID=$(cat /well/nichols/users/scf915/${COHORT}/${COHORT}_subid.txt | awk {'print $1'} | sed "${s}q;d")
	SubID=$(cat /well/nichols/users/scf915/${COHORT}/${COHORT}_subid.txt | awk {'print $1'} | sed "${s}q;d")

        echo $COHORT $TRs $SubID $SesID

        GetAtlas $COHORT $TRs $SubID $SesID

done
