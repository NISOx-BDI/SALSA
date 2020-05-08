#!/bin/bash
set -e

#module load fsl/6.0.3
#module load Python/3.6.6-foss-2018b

module load Python/3.6.6-foss-2018b
module load fsl/6.0.3

COHORT=$1
SubID=$2
SesID="DS2"

#example:
# BIDS_FEAT_prefilt.sh ROCKLAND A00028429 DS2 645

#TR=$4

betflag=1
flag_feat1=1
flag_feat2=1
flag_fnirt1=0
flag_icaaroma=1
#module add fsl

coblar=0

ImageDir="/well/nichols/users/scf915/${COHORT}/raw/sub${SubID}/"

#FUNCIMGNAME=sub-${SubID}_ses-${SesID}_task-rest_acq-${AQLAB}_bold
FUNIMGNAME=rest
FUNCIMG=${ImageDir}/func/${FUNIMGNAME}


ANATIMGNAME=mprage_T1
ANATIMG=${ImageDir}/anat/${ANATIMGNAME}

#cp ${ImageDir}/anat/ ${ImageDir}/anat/${ANATIMGNAME}_brain

STANIMG=${FSLDIR}/data/standard/MNI152_T1_2mm_brain

echo ""
echo "FUNCIMG: ${FUNCIMG}"
echo "ANATIMG: ${ANATIMG}"
echo "STANDARD: ${STANIMG}"

fMRIPREP=/well/nichols/users/scf915/${COHORT}/R_mpp/sub-${SubID}/ses-${SesID}/${FUNIMGNAME}_mpp

if [ $flag_feat1 == 1 ]; then
	[ $coblar == 1 ] && rm -rf ${fMRIPREP}
fi
mkdir -p ${fMRIPREP}

cd ${fMRIPREP}

echo "RESULTS: ${fMRIPREP}"
echo ""

pwd

if [ $betflag == 1 ]; then
	echo ""
	echo "BET the image..."

	${FSLDIR}/bin/robustfov -i ${ImageDir}/anat/mprage_anonymized.nii.gz -r ${ANATIMG}.nii.gz
	${FSLDIR}/bin/bet ${ANATIMG}.nii.gz ${fMRIPREP}/${ANATIMGNAME}_brain.nii.gz -f 0.3 -n -m -R -S -B
fi

if [ $flag_feat1 == 1 ]; then
	echo "Runing the first bit of feat..."

	${FSLDIR}/bin/fslmaths ${FUNCIMG} prefiltered_func_data -odt float
	${FSLDIR}/bin/fslroi prefiltered_func_data example_func 112 1

	${FSLDIR}/bin/mainfeatreg \
	-F 6.00 \
	-d ${fMRIPREP} \
	-l ${fMRIPREP}/feat2_pre \
	-R ${fMRIPREP}/report_unwarp.html \
	-r ${fMRIPREP}/report_reg.html  \
	-i ${fMRIPREP}/example_func.nii.gz  \
	-n 10 \
	-h ${ANATIMGNAME}_brain \
	-w BBR -x 90 \
	-s ${STANIMG} -y 12 -z 90 
else
	echo "NO BBR"
fi

if [ $flag_fnirt1 == 1 ]; then
#runFlirt highres standard $standardDof $standardSearch trilinear $htmlReport "" "$nonLinearResolution" $logFile $fnirtConfig
	    ini=highres
            refi=standard
	    wr=10

	    flirt 

	    immv ${ini}2${refi} ${ini}2${refi}_linear

            ${FSLDIR}/bin/fslmaths $ANATIMG ${ini}_head
	    ${FSLDIR}/bin/fnirt --iout=${ini}2${refi}_head --in=${ini}_head --aff=${ini}2${refi}.mat \
	--cout=${ini}2${refi}_warp --iout=${ini}2${refi} --jout=${ini}2${refi}_jac --config="T1_2_MNI152_2mm" \
	--ref=${refi}_head --refmask=${refi}_mask --warpres=${wr},${wr},${wr}

	    ${FSLDIR}/bin/applywarp -i ${ini} -r ${refi} -o ${ini}2${refi} -w ${ini}2${refi}_warp
else
	echo "Not going to do an independent nonlin registration."
fi


if [ $flag_feat2 == 1 ]; then

	echo "running the second bit of feat..."

	${FSLDIR}/bin/mcflirt -in prefiltered_func_data -out prefiltered_func_data_mcf -mats -plots -reffile example_func -rmsrel -rmsabs -spline_final

	${FSLDIR}/bin/fsl_tsplot -i prefiltered_func_data_mcf.par -t 'MCFLIRT estimated rotations (radians)' -u 1 --start=1 --finish=3 -a x,y,z -w 640 -h 144 -o rot.png
	${FSLDIR}/bin/fsl_tsplot -i prefiltered_func_data_mcf.par -t 'MCFLIRT estimated translations (mm)' -u 1 --start=4 --finish=6 -a x,y,z -w 640 -h 144 -o trans.png
	${FSLDIR}/bin/fsl_tsplot -i prefiltered_func_data_mcf_abs.rms,prefiltered_func_data_mcf_rel.rms -t 'MCFLIRT estimated mean displacement (mm)' -u 1 -w 640 -h 144 -a absolute,relative -o disp.png

	#/bin/mkdir -p mc
	#/bin/mv -f prefiltered_func_data_mcf.mat prefiltered_func_data_mcf.par prefiltered_func_data_mcf_abs.rms prefiltered_func_data_mcf_abs_mean.rms prefiltered_func_data_mcf_rel.rms prefiltered_func_data_mcf_rel_mean.rms mc

	${FSLDIR}/bin/fslmaths prefiltered_func_data_mcf -Tmean mean_func

	${FSLDIR}/bin/bet2 mean_func mask -f 0.3 -n -m
	${FSLDIR}/bin/immv mask_mask mask
	${FSLDIR}/bin/fslmaths prefiltered_func_data_mcf -mas mask prefiltered_func_data_bet

	rm prefiltered_func_data_mcf.nii.gz prefiltered_func_data.nii.gz

	#/usr/local/fsl/bin/fslstats prefiltered_func_data_bet -p 2 -p 98
	#-0.000000 1384.106323
	#${FSLDIR}/bin/fslmaths prefiltered_func_data_bet -thr 138.4106323 -Tmin -bin mask -odt char

	#/usr/local/fsl/bin/fslstats prefiltered_func_data_mcf -k mask -p 50
	#${FSLDIR}/bin/fslmaths mask -dilF mask
	#${FSLDIR}/bin/fslmaths prefiltered_func_data_mcf -mas mask prefiltered_func_data_thresh

else
	echo "NO MCFLIRT"

fi


# RUN ICA-AROMA
if [ $flag_icaaroma == 1 ]; then

	echo "*** ICA-AROMA"

	#module load Python
	rm -rf ${fMRIPREP}/ica-aroma
	mkdir -p ${fMRIPREP}/ica-aroma

	python3.6 /users/nichols/scf915/ICA-AROMA/ICA_AROMA.py \
	-in ${fMRIPREP}/prefiltered_func_data_bet.nii.gz \
	-affmat ${fMRIPREP}/reg/example_func2standard.mat \
	-mc ${fMRIPREP}/prefiltered_func_data_mcf.par \
	-out ${fMRIPREP}/ica-aroma \
	-den both
fi


echo "DONE."
