#bin/bash

set -e

module load fsl

#COHORT=ROCKLAND
#T=240
#TRs=0.645 #2.5 # in second
#NSUB=25
#TaskName=CHECKERBOARD
#TR=$(echo $TRs*1000 | bc | awk -F'.' {'print $1'})
#EDtype="task-${TaskName}_acq-${TR}"

####
#COHORT=tHCP
#T=284
#TRs=0.72 #2.5 # in second
#NSUB=25
#TaskName=MOTOR
#StimulName=lh
#TR=$(echo $TRs*1000 | bc | awk -F'.' {'print $1'})
#EDtype="task-${TaskName}_acq-${TR}_${StimulName}"
######################################
COHORT=tHCP
T=253
TRs=0.72 #2.5 # in second
NSUB=25
TaskName="GAMBLING"
TR=$(echo $TRs*1000 | bc | awk -F'.' {'print $1'})
EDtype="task-${TaskName}_acq-${TR}"
######################################


GSRFLAG=0
ICAFLAG=0

FWHMsize=0

QSUBFLAG=0

TempTreMethodLIST=(hpf hpfk dct poly) #(hpf hpfc hpfk dct spline poly)


TempTrendP=(100 90 3 2)

METHODLIST=(gFAST AR-W ACF gACFadjxACFadj2tT1S0P5)

#METHODLIST=(AR-W)

ARMODE=(1 2 5 10 20 $(echo "sqrt($T)" | bc))

if [ $COHORT == "Cambridge" ];then
	ACMODE=(5 $(echo "sqrt($T)" | bc) $(echo "2*sqrt($T)" | bc) -2)

elif [ $COHORT == "Beijing" ]; then
	ACMODE=(5 10 $(echo "sqrt($T)" | bc) $(echo "2*sqrt($T)" | bc) -2)

elif [ $COHORT == "ROCKLAND" ] && [ $T == 240 ] && [ $TaskName == "CHECKERBOARD" ]; then
	ACMODE=(5 10 $(echo "sqrt($T)" | bc) $(echo "2*sqrt($T)" | bc) -2)
	echo "))))))))))) ROCKLAND TASK ((((((((((((( "
	TempTrendP=(100 90 3 2)

elif [ $COHORT == "tHCP" ] && [ $TaskName == "MOTOR" ]; then
        ACMODE=(5 10 $(echo "sqrt($T)" | bc) $(echo "2*sqrt($T)" | bc) -2)
	TempTrendP=(100 90 4 2)
        echo "))))))))))) MOTOR ((((((((((((( "

elif [ $COHORT == "tHCP" ] && [ $TaskName == "GAMBLING" ]; then
        ACMODE=(5 10 $(echo "sqrt($T)" | bc) $(echo "2*sqrt($T)" | bc) -2)
        TempTrendP=(100 90 3 2)
        echo "))))))))))) GAMBLING ((((((((((((( "
else
	ACMODE=(5 10 15 $(echo "sqrt($T)" | bc) $(echo "2*sqrt($T)" | bc) -2)
fi

if [ $COHORT == "NEO" ]; then
	COHORTDIR=/well/nichols/users/kfh142/data/baby/neofmri_2nd_release_rerun2/
	COHORTPWDIR=/well/nichols/users/scf915/${COHORT}
else
	STRGDIR=/well/nichols/users/scf915
	COHORTDIR=${STRGDIR}/${COHORT}
	COHORTPWDIR=${COHORTDIR}
fi

############################################

if [ $COHORT == "Beijing" ] || [ $COHORT == "ROCKLAND" ]; then
	NCORES=2
elif [ $COHORT == "HCP" ] || [ $COHORT == "NEO" ] || [ $COHORT == "tHCP" ]; then
	NCORES=4
fi

############################################
#/users/nichols/scf915/bin/FILM2/RealfMRI
SUBMITDIR=/users/nichols/scf915/bin/FILM2/RealfMRI/${COHORT}_${T}_${TR}_${TaskName}_gsr${GSRFLAG}_aroma${ICAFLAG}_FLAME1_Submitters
mkdir -p ${SUBMITDIR}

NumJob=0
tr_cnt=0
for TempTreMethod in ${TempTreMethodLIST[@]}
do

	for METH_ID in ${METHODLIST[@]}
	do
		echo ${METH_ID}

		MAO=0; ARMODE0=${ARMODE[@]}

		[ $METH_ID == "ARMAHR" ]&& MAO=1

		if [[ $METH_ID == *"FASTFEAT"* ]] || [[ $METH_ID == *"ACF"* ]]; then
			#ARMODE0=${ACMODE[@]}
			ARMODE0=${ACMODE[@]}
		elif [ $METH_ID == "gFAST" ]; then
			ARMODE0=1
		fi

		for ARO in ${ARMODE0[@]}
		do

			JobName=${COHORT}_${TR}_${T}_${TaskName}_${METH_ID}_AR-${ARO}_MA-${MAO}_FWHM${FWHMsize}_${TempTreMethod}_gsr${GSRFLAG}_aroma${ICAFLAG}
			SubmitterFileName="${SUBMITDIR}/SubmitMe_${JobName}.sh"

			echo ${SubmitterFileName}

			Path2ImgResults=${COHORTPWDIR}/R.PW/${JobName}
			OpLog=${Path2ImgResults}/logs

			mkdir -p ${OpLog}

############################################
############################################
cat > $SubmitterFileName << EOF
#!/bin/bash
#$ -cwd
#$ -q short.qe # short.qc@@short.hge #short.qe    # short.qc@@short.hge
#$ -pe shmem ${NCORES}
#$ -o ${OpLog}/${JobName}_FLAME1_\\\$JOB_ID.out
#$ -e ${OpLog}/${JobName}_FLAME1_\\\$JOB_ID.err
#$ -N ${JobName}

export OMP_NUM_THREADS=${NCORES}

STATFILE=${OpLog}/${JobName}_\${JOB_ID}_\${SGE_TASK_ID}.stat

# The stat file
echo 0 > \$STATFILE

# This whole business is rubbish! This should be fixed!
source \${HOME}/.bashrc
module load fsl

OCTSCRPT=\${HOME}/bin/FILM2/RealfMRI
cd \${OCTSCRPT}

echo \${OCTSCRPT}

##
ConDir="/well/nichols/users/scf915/${COHORT}/R.PW/${COHORT}_${TR}_${T}_${TaskName}_${METH_ID}_AR-${ARO}_MA-0_FWHM${FWHMsize}_${TempTreMethod}_gsr0_aroma0"

LEVEL2DIR=\${ConDir}/LEVEL2

#mkdir -p \${LEVEL2DIR}/FLAME1
rm -r \${LEVEL2DIR}/FLAME1 # clean now, so when the flameo gets to it, doesnt make a duplicate
in4dvol=\${LEVEL2DIR}/ED${EDtype}_${METH_ID}_AR${ARO}_MA0_FWHM${FWHMsize}_${TempTreMethod}${TempTrendP[$tr_cnt]}_Wcbhat_4D
in4dvarcopevol=\${LEVEL2DIR}/ED${EDtype}_${METH_ID}_AR${ARO}_MA0_FWHM${FWHMsize}_${TempTreMethod}${TempTrendP[$tr_cnt]}_wvc_4D
#out4dvol=\${LEVEL2DIR}/FLAME1/ED${EDtype}_${METH_ID}_AR${ARO}_MA0_FWHM${FWHMsize}_${TempTreMethod}${TempTrendP[$tr_cnt]}_Wcbhat
randlog=\${LEVEL2DIR}/FLAME1/ED${EDtype}_${METH_ID}_AR${ARO}_MA0_FWHM${FWHMsize}_${TempTreMethod}${TempTrendP[$tr_cnt]}_Wcbhat_flame1.log

echo "FLAME Log: \$randlog"

sh /users/nichols/scf915/bin/FILM2/mis/run_1sample_flame1.sh \$in4dvol \$in4dvarcopevol >> \$randlog

#flameo --cope=EDtask-GAMBLING_acq-720_ACF_AR30_MA0_FWHM0_dct3_Wcbhat_4D.nii.gz \
#--varcope=EDtask-GAMBLING_acq-720_ACF_AR30_MA0_FWHM0_dct3_wvc_4D.nii.gz \
#--runmode=flame1 \
#--dm=design.mat \
#--tc=cont.mat \
#--mask=mask.nii.gz \
#--cs=design.grp


# The stat file
echo 1 > \$STATFILE

EOF
############################################
############################################
			echo $NumJob
			NumJob=$((NumJob + 1))

			if [ $QSUBFLAG == 1 ]; then
				qsub $SubmitterFileName
			fi
		done
	done
	tr_cnt=$((tr_cnt + 1))
done

echo "Submitted: $NumJob"
