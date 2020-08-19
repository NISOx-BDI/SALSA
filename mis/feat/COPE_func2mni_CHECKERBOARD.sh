#bin/bash

set -e

COHORT=ROCKLAND
T=240
TRs=0.645
TaskName=CHECKERBOARD
NSUB=25

FWHMsize=0

STRGDIR=/well/nichols/users/scf915
COHORTDIR=${STRGDIR}/${COHORT}

TempTreMethodLIST=(dct poly hpf hpfk) #(hpf hpfc hpfk dct spline poly)

TR=$(echo $TRs*1000 | bc | awk -F'.' {'print $1'})

EDtype="task-${TaskName}_acq-${TR}"

#METHODLIST=(ACF gFAST gACFadjxACFadj2tT1S0P5)

METHODLIST=(AR-W)

ARMODE=(1 2 5 10 20)

if [ $COHORT == "ROCKLAND" ] && [ $T == "240" ] && [ $TaskName == "CHECKERBOARD" ]; then
        echo ${COHORT}, ${T}
        TempTrendP=(3 2 100 90)
        ACMODE=(5 10 $(echo "sqrt($T)" | bc) $(echo "2*sqrt($T)" | bc) -2)
else
	echo "Unrecog dataset."
#	exit 1
fi

tr_cnt=0
for TempTreMethod in ${TempTreMethodLIST[@]};do
        echo $tr_cnt
        for METH_ID in ${METHODLIST[@]}; do
                MAO=0; ARMODE0=${ARMODE[@]}
                [ $METH_ID == "ARMAHR" ]&& MAO=1
                if [[ $METH_ID == *"ACF"* ]] ||  [[ $METH_ID == *"FASTFEAT"* ]]; then
                        ARMODE0=${ACMODE[@]}
                elif [ $METH_ID == "gFAST" ]; then
                        ARMODE0=1
                fi

                for ARO in ${ARMODE0[@]}; do

			ConDir="/well/nichols/users/scf915/${COHORT}/R.PW/${COHORT}_${TR}_${T}_${TaskName}_${METH_ID}_AR-${ARO}_MA-0_FWHM${FWHMsize}_${TempTreMethod}_gsr0_aroma0"
                        for (( s=1; s<=${NSUB}; s++ ))
                        do
                                        SubID=$(cat ${COHORTDIR}/${COHORT}_subid.txt | awk {'print $1'} | sed "${s}q;d")
                                        SesID=$(cat ${COHORTDIR}/${COHORT}_sesid.txt | awk {'print $1'} | sed "${s}q;d")


                                        echo "----"
                                        echo ${TempTreMethod}, ${METH_ID}, ${ARO}, ${s}, ${SubID}


                                        SubDIR="${ConDir}/${SubID}_${SesID}"
                                        TDIR=${SubDIR}/ED${EDtype}_${METH_ID}_AR${ARO}_MA0_FWHM${FWHMsize}_${TempTreMethod}${TempTrendP[$tr_cnt]}_cbhat_ICACLEAN0_GSR0.nii.gz
					wTDIR=${SubDIR}/ED${EDtype}_${METH_ID}_AR${ARO}_MA0_FWHM${FWHMsize}_${TempTreMethod}${TempTrendP[$tr_cnt]}_Wcbhat_ICACLEAN0_GSR0.nii.gz

                                        if [ ! -f $TDIR ] || [ ! -f $wTDIR ]; then
                                                echo "$TDIR does not exists."
                                                continue 1
                                        fi

					InDir="/well/nichols/users/scf915/${COHORT}/R_mpp/sub-${SubID}/ses-${SesID}/sub-${SubID}_ses-${SesID}_${EDtype}_bold_mpp"

#					echo $InDir

					sh func2mni_CHECKERBOARD.sh $TDIR $InDir
					sh func2mni_CHECKERBOARD.sh $wTDIR $InDir


                                echo "---"
                        done

			#fslmerge for each scenario
			mkdir -p ${ConDir}/TFCE
			${FSLDIR}/bin/fslmerge -t  ${ConDir}/TFCE/ED${EDtype}_${METH_ID}_AR${ARO}_MA0_FWHM${FWHMsize}_${TempTreMethod}${TempTrendP[$tr_cnt]}_cbhat_4D ${ConDir}/*_*/mni/*_cbhat_ICACLEAN0_GSR0_mni.nii.gz
			echo "4D data is: ${ConDir}/TFCE/ED${EDtype}_${METH_ID}_AR${ARO}_MA0_FWHM${FWHMsize}_${TempTreMethod}${TempTrendP[$tr_cnt]}_cbhat_4D"

			${FSLDIR}/bin/fslmerge -t  ${ConDir}/TFCE/ED${EDtype}_${METH_ID}_AR${ARO}_MA0_FWHM${FWHMsize}_${TempTreMethod}${TempTrendP[$tr_cnt]}_Wcbhat_4D ${ConDir}/*_*/mni/*_Wcbhat_ICACLEAN0_GSR0_mni.nii.gz
                        echo "4D data is: ${ConDir}/TFCE/ED${EDtype}_${METH_ID}_AR${ARO}_MA0_FWHM${FWHMsize}_${TempTreMethod}${TempTrendP[$tr_cnt]}_Wcbhat_4D"


                done
        done
        tr_cnt=$((tr_cnt + 1))
done

