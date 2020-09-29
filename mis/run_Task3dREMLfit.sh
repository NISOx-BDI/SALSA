#bin/bash
set -e

COHORT=$1
TRs=$2
T=$3
SubID=$4
SesID=$5
TaskName=$6
StimulName=$7
FWHMsize=$8
GSRFLAG=$9 #useless but leave it
ICAFLAG=${10} #useless but leave it
Path2ImgResults=${11}
EDtype="ER"


#COHORT=tHCP
#TRs=
#T=
#SubID=
#SesID=
#TaskName=
#StimulName=
#FWHMsize=
#GSRFLAG= #useless but leave it
#ICAFLAG= #useless but leave it
#Path2ImgResults=
#EDtype="ER"

GSRFLAG=0
ICAFLAG=0

###########################################################################
###########################################################################
###########################################################################

TR=$(echo $TRs*1000 | bc | awk -F'.' {'print $1'})

AFNI_SPM_singularity_image=/apps/singularity/afni-r-python3-2020-03-26-v1.sif
AFNI_bin=/opt/afni-latest
PATH2AUX=/users/nichols/scf915/bin/FILM2
#----
prefix=''
if [ $FWHMsize != 0 ]; then prefix=_fwhm${FWHMsize}; fi
if [ $StimulName == 0 ]; then StimulName=""; fi

if [ $COHORT == ROCKLAND ]; then
	COHORTDIR=/well/nichols/users/scf915/${COHORT}
	BOLDDIR=${COHORTDIR}/R_mpp/sub-${SubID}/ses-${SesID}/sub-${SubID}_ses-${SesID}_task-rest_acq-${TR}_bold_mpp/
	BOLDIMG=${BOLDDIR}/prefiltered_func_data_bet${prefix}.nii.gz

elif [ $COHORT == tHCP ]; then
	EDtype=task-${TaskName}_acq-${TR}${StimulName};
	COHORTDIR=/well/nichols/users/scf915/${COHORT}
	Path2ImgRaw=${COHORTDIR}/R_mpp/sub-${SubID}/ses-${SesID}
	BOLDDIR=${Path2ImgRaw}/${SubID}_3T_tfMRI_${SesID}_mpp
	BOLDIMG=${BOLDDIR}/prefiltered_func_data_bet${prefix}.nii.gz
	stim_file=${PATH2AUX}/mis/EVs/${COHORT}/${TaskName}/${EDtype}_event.txt
else
	echo "XXXXXSomething is wrong... "
	exit 1
fi

echo $BOLDIMG

tmpdir=$(mktemp -d)
echo $tmpdir

# Only useful if we had a double task against each other
#stim_file=/users/nichols/scf915/bin/FILM2/mis/EVs/${COHORT}/${COHORT}_sub_${SubID}_T${T}_TR${TR}.txt
#stim_file_E1=$tmpdir/ROCKLAND_sub_${SubID}_T${T}_TR${TR}_E1.txt
#stim_file_E2=$tmpdir/ROCKLAND_sub_${SubID}_T${T}_TR${TR}_E2.txt

#break the file into the E1 and E2
#cat ${stim_file} | awk '{print $1}' > ${stim_file_E1}
#cat ${stim_file} | awk '{print $2}' > ${stim_file_E2}

echo ""
echo "Stimuli file: $stim_file"

AFNIRESULTS=$Path2ImgResults/${SubID}_${SesID}
rm -f ${AFNIRESULTS}/*
mkdir -p ${AFNIRESULTS}
cd ${AFNIRESULTS}

echo ""
echo "We are gonna start in: ${AFNIRESULTS}"

###########################################################################
###########################################################################
###########################################################################

# Later, makes this a bit more elegant...

singularity exec --cleanenv -B ${COHORTDIR} $AFNI_SPM_singularity_image \
$AFNI_bin/3dDeconvolve \
-x1D_stop \
-polort A \
-input $BOLDIMG -automask \
-num_stimts 25 \
-stim_file 1 ${stim_file} -stim_label 1 E1 \
-stim_file 2 $BOLDDIR/mp/${COHORT}_SubID-${SubID}_ses-${SesID}_T${T}_TR${TR}_1.txt -stim_base 2    -stim_label 2  MP1 \
-stim_file 3 $BOLDDIR/mp/${COHORT}_SubID-${SubID}_ses-${SesID}_T${T}_TR${TR}_2.txt -stim_base 3    -stim_label 3  MP2 \
-stim_file 4 $BOLDDIR/mp/${COHORT}_SubID-${SubID}_ses-${SesID}_T${T}_TR${TR}_3.txt -stim_base 4    -stim_label 4  MP3 \
-stim_file 5 $BOLDDIR/mp/${COHORT}_SubID-${SubID}_ses-${SesID}_T${T}_TR${TR}_4.txt -stim_base 5    -stim_label 5  MP4 \
-stim_file 6 $BOLDDIR/mp/${COHORT}_SubID-${SubID}_ses-${SesID}_T${T}_TR${TR}_5.txt -stim_base 6    -stim_label 6  MP5 \
-stim_file 7 $BOLDDIR/mp/${COHORT}_SubID-${SubID}_ses-${SesID}_T${T}_TR${TR}_6.txt -stim_base 7    -stim_label 7  MP6 \
-stim_file 8 $BOLDDIR/mp/${COHORT}_SubID-${SubID}_ses-${SesID}_T${T}_TR${TR}_7.txt -stim_base 8    -stim_label 8  MP7 \
-stim_file 9 $BOLDDIR/mp/${COHORT}_SubID-${SubID}_ses-${SesID}_T${T}_TR${TR}_8.txt -stim_base 9    -stim_label 9 MP8 \
-stim_file 10 $BOLDDIR/mp/${COHORT}_SubID-${SubID}_ses-${SesID}_T${T}_TR${TR}_9.txt -stim_base 10  -stim_label 10 MP9 \
-stim_file 11 $BOLDDIR/mp/${COHORT}_SubID-${SubID}_ses-${SesID}_T${T}_TR${TR}_10.txt -stim_base 11 -stim_label 11 MP10 \
-stim_file 12 $BOLDDIR/mp/${COHORT}_SubID-${SubID}_ses-${SesID}_T${T}_TR${TR}_11.txt -stim_base 12 -stim_label 12 MP11 \
-stim_file 13 $BOLDDIR/mp/${COHORT}_SubID-${SubID}_ses-${SesID}_T${T}_TR${TR}_12.txt -stim_base 13 -stim_label 13 MP12 \
-stim_file 14 $BOLDDIR/mp/${COHORT}_SubID-${SubID}_ses-${SesID}_T${T}_TR${TR}_13.txt -stim_base 14 -stim_label 14 MP13 \
-stim_file 15 $BOLDDIR/mp/${COHORT}_SubID-${SubID}_ses-${SesID}_T${T}_TR${TR}_14.txt -stim_base 15 -stim_label 15 MP14 \
-stim_file 16 $BOLDDIR/mp/${COHORT}_SubID-${SubID}_ses-${SesID}_T${T}_TR${TR}_15.txt -stim_base 16 -stim_label 16 MP15 \
-stim_file 17 $BOLDDIR/mp/${COHORT}_SubID-${SubID}_ses-${SesID}_T${T}_TR${TR}_16.txt -stim_base 17 -stim_label 17 MP16 \
-stim_file 18 $BOLDDIR/mp/${COHORT}_SubID-${SubID}_ses-${SesID}_T${T}_TR${TR}_17.txt -stim_base 18 -stim_label 18 MP17 \
-stim_file 19 $BOLDDIR/mp/${COHORT}_SubID-${SubID}_ses-${SesID}_T${T}_TR${TR}_18.txt -stim_base 19 -stim_label 19 MP18 \
-stim_file 20 $BOLDDIR/mp/${COHORT}_SubID-${SubID}_ses-${SesID}_T${T}_TR${TR}_19.txt -stim_base 20 -stim_label 20 MP19 \
-stim_file 21 $BOLDDIR/mp/${COHORT}_SubID-${SubID}_ses-${SesID}_T${T}_TR${TR}_20.txt -stim_base 21 -stim_label 21 MP20 \
-stim_file 22 $BOLDDIR/mp/${COHORT}_SubID-${SubID}_ses-${SesID}_T${T}_TR${TR}_21.txt -stim_base 22 -stim_label 22 MP21 \
-stim_file 23 $BOLDDIR/mp/${COHORT}_SubID-${SubID}_ses-${SesID}_T${T}_TR${TR}_22.txt -stim_base 23 -stim_label 23 MP22 \
-stim_file 24 $BOLDDIR/mp/${COHORT}_SubID-${SubID}_ses-${SesID}_T${T}_TR${TR}_23.txt -stim_base 24 -stim_label 24 MP23 \
-stim_file 25 $BOLDDIR/mp/${COHORT}_SubID-${SubID}_ses-${SesID}_T${T}_TR${TR}_24.txt -stim_base 25 -stim_label 25 MP24 \
-glt 1 ${PATH2AUX}/mis/EVs/AFNIglts/${COHORT}_${TaskName}_T${T}_TR${TR}_AFNIglt.txt -glt_label 1 'E1'

#/users/nichols/scf915/bin/FILM2/mis/EVs/AFNIglts/HCP_AFNIglt.txt

# Copy a readable design matrix for sanity checks
grep -v '^#' ${AFNIRESULTS}/Decon.xmat.1D > Xmat.txt

###########################################################################
###########################################################################
###########################################################################

echo ""
echo "run 3dREMLfit"

singularity exec --cleanenv -B ${COHORTDIR} $AFNI_SPM_singularity_image \
$AFNI_bin/3dREMLfit \
-matrix ${AFNIRESULTS}/Decon.xmat.1D \
-input $BOLDIMG \
-nobout \
-Rbuck ${AFNIRESULTS}/sub-${SubID}_ses-${SesID}_Decon_REML \
-Rvar ${AFNIRESULTS}/sub-${SubID}_ses-${SesID}_Decon_REMLvar \
-Rwherr ${AFNIRESULTS}/sub-${SubID}_ses-${SesID}_Decon_Rwherr \
-tout \
-Rbeta ${AFNIRESULTS}/sub-${SubID}_ses-${SesID}_Decon_wbeta \
-noFDR \
-verb

###########################################################################
###########################################################################
###########################################################################

# Change all to nifti
echo ""
echo "Reorient everything to nifti and cleanup."

# whitened residuals -------------------------------------
singularity exec --cleanenv -B ${COHORTDIR} $AFNI_SPM_singularity_image \
$AFNI_bin/3dresample -orient RPI -inset ${AFNIRESULTS}/sub-${SubID}_ses-${SesID}_Decon_Rwherr+orig.HEAD -prefix ${AFNIRESULTS}/${EDtype}_sub-${SubID}_ses-${SesID}_Decon_Rwherr.nii.gz

# variance
#singularity exec --cleanenv -B ${COHORTDIR} $AFNI_SPM_singularity_image \
#$AFNI_bin/3dresample -orient RPI -inset ${AFNIRESULTS}/sub-${SubID}_ses-${SesID}_Decon_REMLvar+orig.HEAD -prefix ${AFNIRESULTS}/sub-${SubID}_ses-${SesID}_Decon_REMLvar.nii.gz

# AutoReg Coeff -------------------------------------------
singularity exec --cleanenv -B ${COHORTDIR} $AFNI_SPM_singularity_image \
$AFNI_bin/3dcalc -a ${AFNIRESULTS}/sub-${SubID}_ses-${SesID}_Decon_REMLvar+orig"[a]" -prefix ${AFNIRESULTS}/${EDtype}_sub-${SubID}_ses-${SesID}_Decon_REMLvar_AR.nii.gz -datum float -expr a

# Moving Average Coeff
singularity exec --cleanenv -B ${COHORTDIR} $AFNI_SPM_singularity_image \
$AFNI_bin/3dcalc -a ${AFNIRESULTS}/sub-${SubID}_ses-${SesID}_Decon_REMLvar+orig"[b]" -prefix ${AFNIRESULTS}/${EDtype}_sub-${SubID}_ses-${SesID}_Decon_REMLvar_MA.nii.gz -datum float -expr a

# Lambda
singularity exec --cleanenv -B ${COHORTDIR} $AFNI_SPM_singularity_image \
$AFNI_bin/3dcalc -a ${AFNIRESULTS}/sub-${SubID}_ses-${SesID}_Decon_REMLvar+orig"[lam]" -prefix ${AFNIRESULTS}/${EDtype}_sub-${SubID}_ses-${SesID}_Decon_REMLvar_LAM.nii.gz -datum float -expr a

# Standard Deviation
singularity exec --cleanenv -B ${COHORTDIR} $AFNI_SPM_singularity_image \
$AFNI_bin/3dcalc -a ${AFNIRESULTS}/sub-${SubID}_ses-${SesID}_Decon_REMLvar+orig"[StDev]" -prefix ${AFNIRESULTS}/${EDtype}_sub-${SubID}_ses-${SesID}_Decon_REMLvar_wstd.nii.gz -datum float -expr a

# t-stat ----------------------------------------------------
singularity exec --cleanenv -B ${COHORTDIR} $AFNI_SPM_singularity_image \
$AFNI_bin/3dresample -orient RPI -inset  ${AFNIRESULTS}/sub-${SubID}_ses-${SesID}_Decon_REML+orig.HEAD  -prefix ${AFNIRESULTS}/${EDtype}_sub-${SubID}_ses-${SesID}_Decon_REML_wtv.nii.gz

# whitened betas
singularity exec --cleanenv -B ${COHORTDIR} $AFNI_SPM_singularity_image \
$AFNI_bin/3dresample -orient RPI -inset  ${AFNIRESULTS}/sub-${SubID}_ses-${SesID}_Decon_wbeta+orig.HEAD  -prefix ${AFNIRESULTS}/${EDtype}_sub-${SubID}_ses-${SesID}_Decon_wbeta.nii.gz

# clean HEAD and BRIK files
rm ${AFNIRESULTS}/*.HEAD ${AFNIRESULTS}/*.BRIK
