#bin/bash

set -e

######################################
#COHORT=ROCKLAND
#T=240
#TRs=0.645 #2.5 # in second
#NSUB=25
#TaskName=CHECKERBOARD
#StimulName=0
######################################
COHORT=tHCP
T=284
TRs=0.72 #2.5 # in second
NSUB=25
TaskName=MOTOR
StimulName=lh
######################################

GSRFLAG=0
ICAFLAG=0

FWHMsize=5

QSUBFLAG=1

TempTreMethodLIST=(hpf hpfk dct poly) #(hpf hpfc hpfk dct spline poly)

TR=$(echo $TRs*1000 | bc | awk -F'.' {'print $1'})

METHODLIST=(gFAST AR-W ACF gACFadjT1S0P20 gACFadjT1S0P5 gACFadjxACFadj2tT1S0P5)

ARMODE=(1 2 5 10 20 $(echo "sqrt($T)" | bc))

if [ $COHORT == "Cambridge" ];then
	ACMODE=(5 $(echo "sqrt($T)" | bc) $(echo "2*sqrt($T)" | bc) -2)

elif [ $COHORT == "Beijing" ]; then
	ACMODE=(5 10 $(echo "sqrt($T)" | bc) $(echo "2*sqrt($T)" | bc) -2)

elif [ $COHORT == "ROCKLAND" ] && [ $T==240 ] && [ $TaskName=="CHECKERBOARD" ]; then
	ACMODE=(5 10 $(echo "sqrt($T)" | bc) $(echo "2*sqrt($T)" | bc) -2)
	echo "))))))))))) ROCKLAND TASK ((((((((((((( "
elif [ $COHORT == "tHCP" ] && [ $T==284 ] && [ $TaskName=="MOTOR" ]; then
        ACMODE=(5 10 $(echo "sqrt($T)" | bc) $(echo "2*sqrt($T)" | bc) -2)
        echo "))))))))))) ROCKLAND TASK ((((((((((((( "
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

if [ $COHORT == "Beijing" ] || [ $COHORT == "tHCP" ] || [ $COHORT == "ROCKLAND" ]; then
	NCORES=2
elif [ $COHORT == "HCP" ] || [ $COHORT == "NEO" ] ; then
	NCORES=4
fi

############################################
#/users/nichols/scf915/bin/FILM2/RealfMRI
SUBMITDIR=/users/nichols/scf915/bin/FILM2/RealfMRI/${COHORT}_${T}_${TR}_${TaskName}_gsr${GSRFLAG}_aroma${ICAFLAG}_Submitters
mkdir -p ${SUBMITDIR}

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
#$ -o ${OpLog}/${JobName}_\\\$JOB_ID_\\\$TASK_ID.out
#$ -e ${OpLog}/${JobName}_\\\$JOB_ID_\\\$TASK_ID.err
#$ -N ${JobName}
#$ -t 1-${NSUB} #${NUMJB}

export OMP_NUM_THREADS=${NCORES}

STATFILE=${OpLog}/${JobName}_\${JOB_ID}_\${SGE_TASK_ID}.stat

# The stat file
echo 0 > \$STATFILE

# This whole business is rubbish! This should be fixed!
source \${HOME}/.bashrc
# module use -a /apps/eb/skylake/modules/all
module load Octave/4.4.1-foss-2018b
module load fsl/6.0.3

#SubID=\$(cat ${COHORTDIR}/sub.txt | sed "\${SGE_TASK_ID}q;d" )
SubID=\$(cat ${COHORTPWDIR}/${COHORT}_subid.txt | awk {'print \$1'} | sed "\${SGE_TASK_ID}q;d")
SesID=\$(cat ${COHORTPWDIR}/${COHORT}_sesid.txt | awk {'print \$1'} | sed "\${SGE_TASK_ID}q;d")
#SubID=\$(echo \${SubID} | awk -F"-" '{print \$2}')

OCTSCRPT=\${HOME}/bin/FILM2/RealfMRI
cd \${OCTSCRPT}

echo \${OCTSCRPT}

octave -q --eval "COHORT=\"${COHORT}\"; COHORTDIR=\"${COHORTDIR}\"; Path2ImgResults=\"${Path2ImgResults}\"; TaskName=\"${TaskName}\";  StimulName=\"${StimulName}\"; gsrflag=${GSRFLAG} ; icaclean=${ICAFLAG} ; pwdmethod=\"${METH_ID}\"; lFWHM=${FWHMsize}; TR=${TRs}; Mord=${ARO}; MPparamNum=${MAO}; TempTreMethod=\"${TempTreMethod}\"; SubID=\"\${SubID}\"; SesID=\"\${SesID}\"; RealTask_Img_bmrc; quit"

# The stat file
echo 1 > \$STATFILE

EOF
############################################
############################################
			if [ $QSUBFLAG == 1 ]; then
				qsub $SubmitterFileName
			fi
		done
	done
done
