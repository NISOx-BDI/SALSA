#bin/bash

COHORT=ROCKLAND

METHODLIST=(ACF AR-YW AR-W ARMAHR)
ARMODE=(1 2 5 10 20)
MAO=1

STRGDIR=/well/nichols/users/scf915
COHORTDIR=${STRGDIR}/${COHORT}

SUBMITDIR=/users/nichols/scf915/bin/FILM2/NullRealfMRI/Img/NullImgSim/${COHORT}_Submitters
mkdir -p ${SUBMITDIR}

for METH_ID in ${METHODLIST[@]}
do
	for ARO in ${ARMODE[@]}
	do

	JobName=${COHORT}_${METH_ID}_AR-${ARO}
	SubmitterFileName="${SUBMITDIR}/SubmitMe_${COHORT}_${METH_ID}_AR-${ARO}_MA-${MAO}.sh"

echo ${SubmitterFileName}

Path2ImgRaw=${STRGDIR}/${COHORT}/raw
NUMJB=$(cat ${COHORTDIR}/sub.txt | wc -l )
Path2ImgResults=${COHORTDIR}/R.PW/${METH_ID}_AR-${ARO}
OpLog=${Path2ImgResults}/logs/

mkdir -p ${OpLog}

cat > $SubmitterFileName << EOF

#!/bin/bash
#$ -cwd
#$ -o ${OpLog}/${JobName}_\\\$JOB_ID_\\\$TASK_ID.out
#$ -e ${OpLog}/${JobName}_\\\$JOB_ID_\\\$TASK_ID.err
#$ -N ${JobName}
#$ -t 1-${NUMJB}

module load Octave/4.4.1-foss-2018b

octave -q --eval "try, COHORTDIR=${COHORTDIR}; pwdmethod=${METH_ID}; Mord=${ARO}; subidx=\${SGE_TASK_ID}; ; end; quit"

EOF
	done
done
