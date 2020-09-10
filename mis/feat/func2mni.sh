#bin/bash

set -e

invol=$1
ImgDir=$2
#/well/nichols/users/scf915/ROCKLAND/R.PW/ROCKLAND_NTASK_645_900_gACFadjxACFadj2tT1S0_AR--2_MA-0_FWHM0_hpf_gsr0_aroma0/A00028177_DS2/
#NTASK_SNR_50_EDER_10_gACFadjxACFadj2tT1S0_AR-2_MA0_FWHM0_hpf100_war1_ICACLEAN0_GSR0.nii.gz

involname=$(basename $invol .nii.gz)
involdir=$(dirname $invol)

RegDir=${ImgDir}/reg

#---------------------------------------------------------------------------
warpvol=${ImgDir}/reg/example_func2standard_warp.nii.gz
refvol=${ImgDir}/reg/standard.nii.gz
flirtmat=${ImgDir}/reg/example_func2highres.mat

outvoldir=${involdir}/mni
outvol=${outvoldir}/${involname}_mni.nii.gz
#copy atlas ---------------------------------------------------------------------------
mkdir -p ${outvoldir}

# bring the atlas into the EPI space --------------------------------------------------------------------------
${FSLDIR}/bin/applywarp -i $invol -o $outvol -r $refvol -w $warpvol #--premat=$flirtmat #--interp=nn #--premat=$invflirtmat --interp=nn

echo "Done: $outvol"
