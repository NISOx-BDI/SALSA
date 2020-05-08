#!/bin/bash
#$ -q short.qc
#$ -cwd
#$ -o /users/nichols/scf915/bin/FILM2/NullRealfMRI/Img/feat_sub/logs
#$ -e /users/nichols/scf915/bin/FILM2/NullRealfMRI/Img/feat_sub/logs
#$ -N ROCKLAND_FEAT
#$ -t 4-10

source ${HOME}/.bashrc

module add fsl/6.0.3

COHORTID=ROCKLAND
COHORTDIR=/well/nichols/users/scf915/${COHORTID}
SubID=$(cat ${COHORTDIR}/participants.tsv | awk {'print $1'} | sed "${SGE_TASK_ID}q;d")
SesID=DS2
AQLAB=CAP #1400 #645

echo ""
echo "FEAT IS ON: ${COHORTID} ${SubID} ${SesID} ${TR} "
echo ""

FEATFUNC=/users/nichols/scf915/bin/FILM2/mis/feat
${FEATFUNC}/BIDS_FEAT_prefilt.sh ${COHORTID} ${SubID} ${SesID} ${AQLAB}




