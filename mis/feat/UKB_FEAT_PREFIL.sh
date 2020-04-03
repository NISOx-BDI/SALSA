#!/bin/bash
set -e


flag_feat1=1
flag_feat2=1


SubID=A00029304
SesID=DS2
TR=645

coblar=1

ImageDir=${HOME}/Home/PREW/SimulatefMRI/R/sub-${SubID}/ses-${SesID}
FUNCIMGNAME=sub-${SubID}_ses-${SesID}_task-rest_acq-${TR}_bold

fMRIPREP=${HOME}/Home/PREW/SimulatefMRI/R.mpp/sub-${SubID}/ses-${SesID}/${FUNCIMGNAME}.mpp

if [ $flag_feat1 == 1 ]; then
	[ $coblar == 1 ] && rm -r ${fMRIPREP}
fi

mkdir -p ${fMRIPREP}

FUNCIMG=${ImageDir}/func/${FUNCIMGNAME}.nii.gz
ANATIMG=${ImageDir}/anat/sub-A00029304_ses-DS2_T1w_brain
STANIMG=${FSLDIR}/data/standard/MNI152_T1_2mm_brain

echo ""
echo "FUNCIMG: ${FUNCIMG}"
echo "ANATIMG: ${ANATIMG}"
echo "STANDARD: ${STANIMG}"
echo "RESULTS: ${fMRIPREP}"
echo ""

cd ${fMRIPREP}

${FSLDIR}/bin/fslmaths ${FUNCIMG} prefiltered_func_data -odt float
${FSLDIR}/bin/fslroi prefiltered_func_data example_func 112 1

if [ $flag_feat1 == 1 ]; then
	echo "Runing the first bit of feat..."

	${FSLDIR}/bin/mainfeatreg \
	-F 6.00 \
	-d ${fMRIPREP} \
	-l ${fMRIPREP}/feat2_pre \
	-R ${fMRIPREP}/report_unwarp.html \
	-r ${fMRIPREP}/report_reg.html  \
	-i ${fMRIPREP}/example_func.nii.gz  \
	-h ${ANATIMG} \
	-w BBR -x 90 \
	-s ${STANIMG} -y 12 -z 90 
else
	echo "Pass the first stage!"
fi

if [ $flag_feat2 == 1 ]; then 
	
	echo "running the second bit of feat..."

	/usr/local/fsl/bin/mcflirt -in prefiltered_func_data -out prefiltered_func_data_mcf -mats -plots -reffile example_func -rmsrel -rmsabs -spline_final
	
	${FSLDIR}/bin/fsl_tsplot -i prefiltered_func_data_mcf.par -t 'MCFLIRT estimated rotations (radians)' -u 1 --start=1 --finish=3 -a x,y,z -w 640 -h 144 -o rot.png 
	${FSLDIR}/bin/fsl_tsplot -i prefiltered_func_data_mcf.par -t 'MCFLIRT estimated translations (mm)' -u 1 --start=4 --finish=6 -a x,y,z -w 640 -h 144 -o trans.png 
	${FSLDIR}/bin/fsl_tsplot -i prefiltered_func_data_mcf_abs.rms,prefiltered_func_data_mcf_rel.rms -t 'MCFLIRT estimated mean displacement (mm)' -u 1 -w 640 -h 144 -a absolute,relative -o disp.png 

	#/bin/mkdir -p mc  
	#/bin/mv -f prefiltered_func_data_mcf.mat prefiltered_func_data_mcf.par prefiltered_func_data_mcf_abs.rms prefiltered_func_data_mcf_abs_mean.rms prefiltered_func_data_mcf_rel.rms prefiltered_func_data_mcf_rel_mean.rms mc

	${FSLDIR}/bin/fslmaths prefiltered_func_data_mcf -Tmean mean_func

	${FSLDIR}/bin/bet2 mean_func mask -f 0.3 -n -m 
	#${FSLDIR}/bin/immv mask_mask mask

	${FSLDIR}/bin/fslmaths prefiltered_func_data_mcf -mas mask prefiltered_func_data_bet

	#/usr/local/fsl/bin/fslstats prefiltered_func_data_bet -p 2 -p 98
	#-0.000000 1384.106323 
	#${FSLDIR}/bin/fslmaths prefiltered_func_data_bet -thr 138.4106323 -Tmin -bin mask -odt char

	#/usr/local/fsl/bin/fslstats prefiltered_func_data_mcf -k mask -p 50
	#${FSLDIR}/bin/fslmaths mask -dilF mask
	#${FSLDIR}/bin/fslmaths prefiltered_func_data_mcf -mas mask prefiltered_func_data_thresh
fi
