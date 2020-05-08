#!/bin/bash
#$ -q short.qe # short.qc@@short.hge
#$ -cwd
#$ -pe shmem 4
#$ -o /users/nichols/scf915/bin/FILM2/NullRealfMRI/Img/feat_sub/logs
#$ -e /users/nichols/scf915/bin/FILM2/NullRealfMRI/Img/feat_sub/logs
#$ -N FCON1000_FEAT
#$ -t 1-82

source ${HOME}/.bashrc
export OMP_NUM_THREADS=4

COHORTID="Beijing" #Cambridge
COHORTDIR=/well/nichols/users/scf915/${COHORTID}
SubID=$(cat ${COHORTDIR}/${COHORTID}_subid.txt | awk {'print $1'} | sed "${SGE_TASK_ID}q;d" | cut -d'b' -f 2)

echo ""
echo "FEAT IS ON: ${COHORTID} ${SubID}"
echo ""

FEATFUNC=/users/nichols/scf915/bin/FILM2/mis/feat
sh ${FEATFUNC}/FEAT_fcon1000_prefilt.sh ${COHORTID} ${SubID}
# -- -- -- ${FEATFUNC}/BIDS_FEAT_prefilt.sh ${COHORTID} ${SubID} ${SesID} ${AQLAB}




