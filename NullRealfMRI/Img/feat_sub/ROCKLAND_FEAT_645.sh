#!/bin/bash
#$ -q short.qe # short.qc@@short.hge
#$ -cwd
#$ -pe shmem 4
#$ -o /users/nichols/scf915/bin/FILM2/NullRealfMRI/Img/feat_sub/logs
#$ -e /users/nichols/scf915/bin/FILM2/NullRealfMRI/Img/feat_sub/logs
#$ -N ROCKLAND_FEAT
#$ -t 1-60

source ${HOME}/.bashrc
export OMP_NUM_THREADS=4

#module add fsl/6.0.3

COHORTID=ROCKLAND
COHORTDIR=/well/nichols/users/scf915/${COHORTID}
SubID=$(cat ${COHORTDIR}/ROCKLAND_subid.txt | awk {'print $1'} | sed "${SGE_TASK_ID}q;d")
SesID=DS2
AQLAB="task-CHECKERBOARD_acq-645"

#AQLAB=645 #1400 #CAP #645

echo ""
echo "FEAT IS ON: ${COHORTID} ${SubID} ${SesID} ${TR} "
echo ""

FEATFUNC=/users/nichols/scf915/bin/FILM2/mis/feat
${FEATFUNC}/BIDS_FEAT_prefilt.sh ${COHORTID} ${SubID} ${SesID} ${AQLAB}




