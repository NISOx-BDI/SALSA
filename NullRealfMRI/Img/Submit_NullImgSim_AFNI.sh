#bin/bash

set -e

# COHORT      T    TR
# Cambridge  119   3
# Beijing    225   2
# ROCKLAND   900   0.645
# ROCKLAND   120   2.5
# ROCKLAND   404   1.4
# HCP        1200  0.72
# NEO	     2300  0.392

#COHORT=ROCKLAND
#T=404 #120
#TRs=1.4 #2.5 # in second
#NSUB=51
#EDtype="ER"

#COHORT=Beijing
#T=225 #120
#TRs=2 #2.5 # in second
#NSUB=51
#EDtype="ER"

#COHORT=ROCKLAND
#T=900 #120
#TRs=0.645 #2.5 # in second
#NSUB=51
#EDtype="ER"

COHORT=NEO
T=2300 #120
TRs=0.392 #2.5 # in second
NSUB=51
EDtype="ER"

#COHORT=HCP
#T=1200 #120
#TRs=0.72 #2.5 # in second
#NSUB=51
#EDtype="ER"

######################################

GSRFLAG=0
ICAFLAG=0

FWHMsize=5
QSUBFLAG=1

ARO=1
MAO=1
TempTreMethod=poly
METH_ID=3dREMLfit

TR=$(echo $TRs*1000 | bc | awk -F'.' {'print $1'})
############################################
if [ $COHORT == "Beijing" ] || [ $COHORT == "ROCKLAND" ]; then
	NCORES=2
elif [ $COHORT == "HCP" ] || [ $COHORT == "NEO" ] ; then
	NCORES=4
fi
############################################
STRGDIR=/well/nichols/users/scf915
COHORTDIR=${STRGDIR}/${COHORT}
COHORTPWDIR=${COHORTDIR}
SUBMITDIR=/users/nichols/scf915/bin/FILM2/NullRealfMRI/Img/${COHORT}_${T}_${TR}_gsr${GSRFLAG}_aroma${ICAFLAG}_Submitters
mkdir -p ${SUBMITDIR}
JobName=${COHORT}_${TR}_${T}_${METH_ID}_AR-${ARO}_MA-${MAO}_FWHM${FWHMsize}_${TempTreMethod}_gsr${GSRFLAG}_aroma${ICAFLAG}
SubmitterFileName="${SUBMITDIR}/SubmitMe_${JobName}.sh"
echo $SubmitterFileName
Path2ImgResults=${COHORTPWDIR}/R.PW/${JobName}
OpLog=${Path2ImgResults}/logs

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

# Run 3dREMLfit -----
#COHORT=$1 TRs=$2 T=$3 SubID=$4 SesID=$5 FWHMsize=$6 GSRFLAG=$7 ICAFLAG=$8
sh /users/nichols/scf915/bin/FILM2/mis/run_Null3dREMLfit.sh $COHORT $TRs $T \$SubID \$SesID $FWHMsize $GSRFLAG $ICAFLAG $Path2ImgResults

# The stat file
echo 1 > \$STATFILE

EOF

############################################
############################################

if [ $QSUBFLAG == 1 ]; then
	qsub $SubmitterFileName
fi

