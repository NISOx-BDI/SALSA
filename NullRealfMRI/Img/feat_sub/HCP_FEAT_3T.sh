#!/bin/bash
#$ -q short.qe # short.qc@@short.hge
#$ -cwd
#$ -pe shmem 4
#$ -o /users/nichols/scf915/bin/FILM2/NullRealfMRI/Img/feat_sub/logs
#$ -e /users/nichols/scf915/bin/FILM2/NullRealfMRI/Img/feat_sub/logs
#$ -N HCP_FEAT
#$ -t 1-25

source ${HOME}/.bashrc
export OMP_NUM_THREADS=4

#module add fsl/6.0.3

COHORTID=tHCP
COHORTDIR=/well/nichols/users/scf915/${COHORTID}
SubID=$(cat ${COHORTDIR}/${COHORTID}_subid.txt | awk {'print $1'} | sed "${SGE_TASK_ID}q;d")
SesID=GAMBLING_LR

echo ""
echo "FEAT IS ON: ${COHORTID} ${SubID} ${SesID} ${TR} "
echo ""

FEATFUNC=/users/nichols/scf915/bin/FILM2/mis/feat
${FEATFUNC}/HCP_FEAT_prefilt.sh ${COHORTID} ${SubID} ${SesID}




