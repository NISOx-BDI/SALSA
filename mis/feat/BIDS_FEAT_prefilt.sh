#!/bin/bash
set -e

#module load fsl/6.0.3
#module load Python/3.6.6-foss-2018b

module load Python/3.6.6-foss-2018b
module load fsl/6.0.3

COHORT=$1
SubID=$2
SesID=$3
AQLAB=$4

#example:
# BIDS_FEAT_prefilt.sh ROCKLAND A00028429 DS2 645

#TR=$4

betflag=0
flag_feat1=0
flag_fast=1
flag_feat2=0
flag_fnirt1=1
flag_smooth=0
flag_icaaroma=0

#module add fsl

#smoothing parameters
FWHMl=5

coblar=0

ImageDir="/well/nichols/users/scf915/${COHORT}/raw/sub-${SubID}/ses-${SesID}"

#FUNCIMGNAME=sub-${SubID}_ses-${SesID}_task-rest_acq-${AQLAB}_bold
FUNIMGNAME=sub-${SubID}_ses-${SesID}_task-rest_acq-${AQLAB}_bold
FUNCIMG=${ImageDir}/func/${FUNIMGNAME}

ANATIMGNAME=sub-${SubID}_ses-${SesID}_T1w
ANATIMG=${ImageDir}/anat/${ANATIMGNAME}

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

cp ${ANATIMG}.nii.gz ${fMRIPREP}/ # take the image to the R_mpp directory

if [ $betflag == 1 ]; then
	echo ""
	echo "BET the image..."
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
	sh ${HOME}/bin/FILM2/mis/feat/reg_fun2mni.sh ${fMRIPREP}/${ANATIMGNAME}
else
	echo "Not going to do an independent nonlin registration."
fi


if [ $flag_fast == 1 ]; then
	mkdir -p ${fMRIPREP}/seg
	echo "Applying fast on T1 image."

	# Apply FAST on the T1 images
	${FSLDIR}/bin/fast -n 3 -g \
	-o ${fMRIPREP}/seg/${ANATIMGNAME} \
	${fMRIPREP}/${ANATIMGNAME}_brain.nii.gz

       ${FSLDIR}/bin/applywarp --interp=nn \
       --ref=example_func \
       --in=${fMRIPREP}/seg/${ANATIMGNAME}_seg \
       --out=${fMRIPREP}/seg/${ANATIMGNAME}_func_seg \
       --premat=${fMRIPREP}/reg/highres2example_func.mat


	# Take FAST results into the EPI space.
	for segi in 0 1 2; do
		echo "Take ${fMRIPREP}/seg/${ANATIMGNAME}_seg_${segi} into EPI space"

		${FSLDIR}/bin/applywarp --interp=nn \
		--ref=example_func \
		--in=${fMRIPREP}/seg/${ANATIMGNAME}_seg_${segi} \
		--out=${fMRIPREP}/seg/${ANATIMGNAME}_func_seg_${segi} \
		--premat=${fMRIPREP}/reg/highres2example_func.mat

		echo "Take ${fMRIPREP}/seg/${ANATIMGNAME}_pve_${segi} into EPI space"
		${FSLDIR}/bin/applywarp --interp=trilinear \
                --ref=example_func \
                --in=${fMRIPREP}/seg/${ANATIMGNAME}_pve_${segi} \
                --out=${fMRIPREP}/seg/${ANATIMGNAME}_func_pve_${segi} \
                --premat=${fMRIPREP}/reg/highres2example_func.mat

	done
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

# Do the smoothing
if [ $flag_smooth == 1 ]; then

	origimage=${fMRIPREP}/prefiltered_func_data_bet
	finalresult=${origimage}_fwhm${FWHMl}

	maskimage=${fMRIPREP}/mask

	smoothparSig=$(bc -l <<< "${FWHMl}/2.3548")

	echo "===SMOOTHING:"
	echo "-IN:   ${origimage}"
	echo "-MASK: ${maskimage}"
	echo "-OUT:  ${finalresult}"
	echo "-KERNEL(mm): ${FWHMl} , sigma: ${smoothparSig}"
	echo ""


	${FSLDIR}/bin/fslmaths $origimage -s $smoothparSig -mas $maskimage tmp_result1_tmp
	${FSLDIR}/bin/fslmaths $maskimage -s $smoothparSig -mas $maskimage $finalresult
	${FSLDIR}/bin/fslmaths tmp_result1_tmp -div $finalresult $finalresult

	PREFILTIMG=${finalresult}

	${FSLDIR}/bin/imrm tmp_result1_tmp
else
	PREFILTIMG=${fMRIPREP}/prefiltered_func_data_bet
	echo "NO SMOOTHING"
fi



# RUN ICA-AROMA
if [ $flag_icaaroma == 1 ]; then

	echo "*** ICA-AROMA"

	#module load Python
	ica_dir=${fMRIPREP}/ica-aroma_fwhm${FWHMl}
	rm -rf ${ica_dir}
	mkdir -p ${ica_dir}

	python3.6 /users/nichols/scf915/ICA-AROMA/ICA_AROMA.py \
	-in ${PREFILTIMG}.nii.gz \
	-affmat ${fMRIPREP}/reg/example_func2standard.mat \
	-mc ${fMRIPREP}/prefiltered_func_data_mcf.par \
	-out ${ica_dir} \
	-den both
fi


echo "DONE."
