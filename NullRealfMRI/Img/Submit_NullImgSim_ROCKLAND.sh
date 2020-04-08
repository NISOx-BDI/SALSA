#bin/bash

COHORT=ROCKLAND

METHODLIST=(ACF AR-YW AR-W ARMAHR)
ARMODE=(1 2 5 10 20)
FWHMsize=0

STRGDIR=/well/nichols/users/scf915
COHORTDIR=${STRGDIR}/${COHORT}

Path2ImgRaw=${STRGDIR}/${COHORT}/raw
TR=$(cat ${Path2ImgRaw}/task-rest_acq-645_bold.json | grep RepetitionTime | awk {'print $2'} | awk -F"," '{print $1}')

SesID=DS2

NUMJB=$(cat ${COHORTDIR}/sub.txt | wc -l )

SUBMITDIR=/users/nichols/scf915/bin/FILM2/NullRealfMRI/Img/${COHORT}_Submitters
mkdir -p ${SUBMITDIR}

for METH_ID in ${METHODLIST[@]}
do
	MAO=0
	[ $METH_ID == "ARMAHR" ]&& MAO=1

	for ARO in ${ARMODE[@]}
	do

	JobName=${COHORT}_${METH_ID}_AR-${ARO}
	SubmitterFileName="${SUBMITDIR}/SubmitMe_${COHORT}_${METH_ID}_AR-${ARO}_MA-${MAO}_FWHM${FWHMsize}.sh"

echo ${SubmitterFileName}

Path2ImgResults=${COHORTDIR}/R.PW/${METH_ID}_AR-${ARO}_FWHM${FWHMsize}
OpLog=${Path2ImgResults}/logs/


mkdir -p ${OpLog}

cat > $SubmitterFileName << EOF
#!/bin/bash
#$ -cwd
#$ -q short.qc@@short.hge
#$ -o ${OpLog}/${JobName}_\\\$JOB_ID_\\\$TASK_ID.out
#$ -e ${OpLog}/${JobName}_\\\$JOB_ID_\\\$TASK_ID.err
#$ -N ${JobName}
#$ -t 1-200 #${NUMJB}

# This whole business is rubbish! This should be fixed!
# source \${HOME}/.bashrc
# module use -a /apps/eb/skylake/modules/all
module load Octave/4.4.1-foss-2018b

SubID=\$(cat ${COHORTDIR}/sub.txt | sed "\${SGE_TASK_ID}q;d" )
SubID=\$(echo \${SubID} | awk -F"-" '{print \$2}')

OCTSCRPT=\${HOME}/bin/FILM2/NullRealfMRI/Img
cd \${OCTSCRPT}
octave -q --eval "COHORTDIR=\"${COHORTDIR}\"; pwdmethod=\"${METH_ID}\"; lFWHM=${FWHMsize}; TR=${TR}; Mord=${ARO}; SubID=\"\${SubID}\"; SesID=\"${SesID}\"; NullSim_Img_bmrc; quit"

EOF
	done
done
