#bin/bash


#function ToStandard () {

DataDir=/well/nichols/users/scf915

invol=$1
#/well/nichols/users/scf915/ROCKLAND/R.PW/ROCKLAND_NTASK_645_900_gACFadjxACFadj2tT1S0_AR--2_MA-0_FWHM0_hpf_gsr0_aroma0/A00028177_DS2/
#NTASK_SNR_50_EDER_10_gACFadjxACFadj2tT1S0_AR-2_MA0_FWHM0_hpf100_war1_ICACLEAN0_GSR0.nii.gz

involname=$(basename $invol .nii.gz)
involdir=$(dirname $invol)

# Get info out of the input

Cohort=$(echo $involdir | awk -F"/" '{print $6}')

SubIDSesID=$(echo $involdir | awk -F"/" '{print $9}')
SubID=$(echo $SubIDSesID | awk -F"_" '{print $1}')
SesID=$(echo $SubIDSesID | awk -F"_" '{print $2}')

DIRNAME=$(echo $involdir | awk -F"/" '{print $8}')
TR=$(echo $DIRNAME | awk -F"_" '{print $3}')

SNR=$(echo $involname | awk -F"_" '{print $3}')
ER=$(echo $involname | awk -F"_" '{print $4}')
BCl=$(echo $involname | awk -F"_" '{print $5}')
PWM=$(echo $involname | awk -F"_" '{print $6}')
ARM=$(echo $involname | awk -F"_" '{print $7}')
MAM=$(echo $involname | awk -F"_" '{print $8}')
FWHMl=$(echo $involname | awk -F"_" '{print $9}')
DTR=$(echo $involname | awk -F"_" '{print $10}')
IMT=$(echo $involname | awk -F"_" '{print $11}')
ICA=$(echo $involname | awk -F"_" '{print $12}')
GSR=$(echo $involname | awk -F"_" '{print $13}')

#


Path2ImgDir=${DataDir}/${Cohort}/R_mpp

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
refvol=${ImgDir}/reg/standard.nii.gz
flirtmat=${ImgDir}/reg/example_func2highres.mat

outvoldir=${involdir}/mni
outvol=${outvoldir}/NTASK_SNR_${SNR}_${ER}_${BCl}_${PWM}_${ARM}_${MAM}_${FWHMl}_${DTR}_${IMT}_${ICA}_${GSR}_mni.nii.gz
#copy atlas ---------------------------------------------------------------------------
mkdir -p ${outvoldir}

# bring the atlas into the EPI space --------------------------------------------------------------------------
${FSLDIR}/bin/applywarp -i $invol -o $outvol -r $refvol -w $warpvol #--premat=$flirtmat #--interp=nn #--premat=$invflirtmat --interp=nn

#}



#MAIN ######################################################################
############################################################################
#
#invol=$1
#
#TR=$(echo $TRs*1000 | bc | awk -F'.' {'print $1'})
#
#for (( s=1; s<=51; s++ ))
#do
#        #SubID=$(cat /well/nichols/users/scf915/${COHORT}/${COHORT}_subid.txt | awk {'print $1'} | sed "${s}q;d")
#	SubID=$(cat /well/nichols/users/scf915/${COHORT}/${COHORT}_subid.txt | awk {'print $1'} | sed "${s}q;d")
#
#        echo $COHORT $TRs $SubID $SesID
#        GetAtlas $invol
#
#done
