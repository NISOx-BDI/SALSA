#bin/bash

set -e

NCORES=2

# COHORT      T    TR
# Cambridge  119   3
# Beijing    225   2
# ROCKLAND   900   0.645
# ROCKLAND   120   2.5
# ROCKLAND   404   1.4
# HCP        1200  0.72

#COHORT=ROCKLAND
#T=900 #120
#TRs=0.645 #2.5 # in second

COHORT=NEO
T=2300
TRs=0.392
NUMJB=51

GSRFLAG=0
ICAFLAG=0

FWHMsize=0

EDtype="ERF"

QSUBFLAG=1

TempTreMethodLIST=(hpf hpfk dct poly) #(hpf hpfc hpfk dct spline poly)

TR=$(echo $TRs*1000 | bc | awk -F'.' {'print $1'})

#ACF ACFadj AR-W ARMAHR
#METHODLIST=(ACF ACFadj gFASTxACFadj2t2j) #ACFT0S0 ACFadjT1S0 ACFadjT0S5 ACFadjT0S0 gFAST gFASTxACFadj2t2jT0S0) #(FASTFEAT20T0S0 FASTFEAT50T0S0 FASTFEAT100T0S0) #(ACFT0S0 ACFadjT1S0 ACFadjT0S5 ACFadjT0S0 gFASTxACFadj2t2jT0S0) #(ACFadjT1S0) #(ACFadjT0S5 ACFT0S0 FASTFEAT20T0S5 FASTFEAT50T0S5 FASTFEAT100T0S5)  #(gFASTxACFadj2t2jT0S0) #(ACFadjT0S0 gFASTxACFadj2t2jT0S0 FASTFEAT20T0S0 FASTFEAT50T0S0 FASTFEAT100T0S0)  #(FASTFEAT10 FASTFEAT20 FASTFEAT50 FASTFEAT100)  #(vFAST20 vFAST50 vFAST100) #AR-W #gFASTxACFadj #(ACF ACFadj gFASTxACFadj2t gFASTxACFadj2t2j AR-W) #(gReMLxARw2t2j) #(ACF ACFadj gFASTxACFadj2t AR-W) #gFASTxACFadj2t2j #gFASTxACFadj2t #gFAST #(ACF ACFadj AR-W) #(ACF ACFadj ACFadjx2 AR-W) #ARMAHR)

# The global ACFadj family of methods
#METHODLIST=(gACFadjT1S0 gACFadjxACFadjT1S0 gACFadjxACFadj2tT1S0)

# Methods for Task
#METHODLIST=(gACFadjT1S0 gFAST ACFadj ACF gACFadjxACFadj2tT1S0 AR-W)
#METHODLIST=(gACFadjxACFadj2tT1S0)
METHODLIST=(gFAST AR-W ACF gACFadjxACFadj2tT1S0P5)

ARMODE=(1 2 5 10 20 $(echo "sqrt($T)" | bc))

if [ $COHORT == "Cambridge" ];then
	ACMODE=(5 $(echo "sqrt($T)" | bc) $(echo "2*sqrt($T)" | bc) -2)

elif [ $COHORT == "Beijing" ]; then
	ACMODE=(5 10 $(echo "sqrt($T)" | bc) $(echo "2*sqrt($T)" | bc) -2)
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
#TR=$(cat ${Path2ImgRaw}/task-rest_acq-645_bold.json | grep RepetitionTime | awk {'print $2'} | awk -F"," '{print $1}')

#if [ $COHORT == "Beijing" ] || [ $COHORT == "ROCKLAND" ]; then
#	SesID=DS2
#elif [ $COHORT == "HCP" ]; then
#	SesID=REST1_LR
#fi

#NUMJB=$(cat ${COHORTDIR}/sub.txt | wc -l )
############################################

SUBMITDIR=/users/nichols/scf915/bin/FILM2/NullRealfMRI/Img/${COHORT}_NTASK_${T}_${TR}_gsr${GSRFLAG}_aroma${ICAFLAG}_Submitters
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

			JobName=${COHORT}_NTASK_${TR}_${T}_${METH_ID}_AR-${ARO}_MA-${MAO}_FWHM${FWHMsize}_${TempTreMethod}_gsr${GSRFLAG}_aroma${ICAFLAG}
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
#$ -t 1-${NUMJB}

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

OCTSCRPT=\${HOME}/bin/FILM2/NullRealfMRI/Img
cd \${OCTSCRPT}
octave -q --eval "COHORT=\"${COHORT}\" ;COHORTDIR=\"${COHORTDIR}\"; Path2ImgResults=\"${Path2ImgResults}\"; EDtype=\"${EDtype}\"; gsrflag=${GSRFLAG} ; icaclean=${ICAFLAG} ; pwdmethod=\"${METH_ID}\"; lFWHM=${FWHMsize}; TR=${TRs}; Mord=${ARO}; MPparamNum=${MAO}; TempTreMethod=\"${TempTreMethod}\"; SubID=\"\${SubID}\"; SesID=\"\${SesID}\"; NullTask_Img00_bmrc; quit"

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
