#!/bin/bash
#$ -q short.qe # short.qc@@short.hge
#$ -cwd
#$ -pe shmem 4
#$ -o /users/nichols/scf915/bin/FILM2/NullRealfMRI/Img/feat_sub/logs
#$ -e /users/nichols/scf915/bin/FILM2/NullRealfMRI/Img/feat_sub/logs
#$ -N NEO_FEAT
#$ -t 1-60

source ${HOME}/.bashrc
export OMP_NUM_THREADS=4

module add fsl/6.0.3

COHORT=NEO

SubID=$(cat /well/nichols/users/scf915/${COHORT}/${COHORT}_subid.txt | awk {'print $1'} | sed "${SGE_TASK_ID}q;d")
SesID=$(cat /well/nichols/users/scf915/${COHORT}/${COHORT}_sesid.txt | awk {'print $1'} | sed "${SGE_TASK_ID}q;d")
echo $COHORT $SubID $SesID
sh /users/nichols/scf915/bin/FILM2/mis/feat/neofmri_prep.sh $SubID $SesID


