#!bin/bash

#sub-A00028185_ses-DS2
SubID=A00028185
SesID=DS2

INPUTDIR=/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/ROCKLAND/sub-${SubID}/ses-${SesID}/sub-${SubID}_ses-${SesID}_task-rest_acq-645_bold_mpp

cd ${INPUTDIR}
INPUTIMG=${INPUTDIR}/sub-${SubID}_ses-${SesID}_T1w_brain.nii.gz

RESULTDIR=${INPUTDIR}/seg
mkdir -p ${RESULTDIR}

RESULTIMG=${RESULTDIR}/sub-${SubID}_ses-${SesID}

${FSLDIR}/bin/fast -n 3 \
--iter=10 \
-f 0.1 \
--segments -B -p \
--out=${RESULTIMG} ${INPUTIMG}

${FSLDIR}/bin/flirt -in seg/sub-A00028185_ses-DS2_pve_1.nii.gz \
-ref mean_func.nii.gz \
-applyxfm \
-interp nearestneighbour \
-init epi2struct/struct2epi.mat \
-out GMepi2

${FSLDIR}/bin/fslmaths GMepi2 -thr .20 -bin GMepi2mask

${FSLDIR}/bin/fslmaths prefiltered_func_data_bet.nii.gz -mas GMepi2mask prefiltered_func_data_gm2_bet.nii.gz
