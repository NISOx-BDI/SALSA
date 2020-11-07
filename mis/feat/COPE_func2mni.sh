#bin/bash

set -e

COHORT=ROCKLAND
T=240
TRs=0.645
TaskName=CHECKERBOARD
NSUB=25
StimulName=""
####################################
#COHORT=tHCP
#T=284
#TRs=0.72
#TaskName=MOTOR
#NSUB=25
#StimulName="_lh"
###################################
#COHORT=tHCP
#T=253
#TRs=0.72
#TaskName=GAMBLING
#NSUB=25
#StimulName=""
###################################

FWHMsize=5

STRGDIR=/well/nichols/users/scf915
COHORTDIR=${STRGDIR}/${COHORT}

TempTreMethodLIST=(dct poly hpf hpfk) #(hpf hpfc hpfk dct spline poly)

TR=$(echo $TRs*1000 | bc | awk -F'.' {'print $1'})

EDtype="task-${TaskName}_acq-${TR}${StimulName}"
EDtype_Rmpp="task-${TaskName}_acq-${TR}"

#METHODLIST=(ACF gFAST gACFadjxACFadj2tT1S0P5 AR-W)

METHODLIST=(3dREMLfit)

ARMODE=(1 2 5 10 20)

if [ $COHORT == "ROCKLAND" ] && [ $T == "240" ] && [ $TaskName == "CHECKERBOARD" ]; then
        echo ${COHORT}, ${T}
        TempTrendP=(3 2 100 90)
        ACMODE=(5 10 $(echo "sqrt($T)" | bc) $(echo "2*sqrt($T)" | bc) -2)
elif [ $COHORT == "tHCP" ] && [ $TaskName == "MOTOR" ]; then
	TempTrendP=(4 2 100 90)
        ACMODE=(5 10 $(echo "sqrt($T)" | bc) $(echo "2*sqrt($T)" | bc) -2)
elif [ $COHORT == "tHCP" ] && [ $TaskName == "GAMBLING" ]; then
	TempTrendP=(3 2 100 90)
        ACMODE=(5 10 $(echo "sqrt($T)" | bc) $(echo "2*sqrt($T)" | bc) -2)
else
	echo "Unrecog dataset."
#	exit 1
fi

if [ ${#METHODLIST[@]} == 1 ] & [ $METHODLIST == '3dREMLfit' ]; then
	TempTreMethodLIST=(poly)
	TempTrendP=2
#	if [ $COHORT == "ROCKLAND" ] && [ $T == "240" ] && [ $TaskName == "CHECKERBOARD" ]; then
#        	echo ${COHORT}, ${T}
#        	TempTrendP=(2)
#	elif [ $COHORT == "tHCP" ] && [ $TaskName == "MOTOR" ]; then
#        	TempTrendP=(2)
#	elif [ $COHORT == "tHCP" ] && [ $TaskName == "GAMBLING" ]; then
#        	TempTrendP=(2)
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
		elif [ $METH_ID == "3dREMLfit" ]; then
			ARMODE0=1
			MAO=1
                fi

                for ARO in ${ARMODE0[@]}; do

			ConDir="/well/nichols/users/scf915/${COHORT}/R.PW/${COHORT}_${TR}_${T}_${TaskName}_${METH_ID}_AR-${ARO}_MA-${MAO}_FWHM${FWHMsize}_${TempTreMethod}_gsr0_aroma0"
                        for (( s=1; s<=${NSUB}; s++ ))
                        do
                                        SubID=$(cat ${COHORTDIR}/${COHORT}_subid.txt | awk {'print $1'} | sed "${s}q;d")
                                        SesID=$(cat ${COHORTDIR}/${COHORT}_${TaskName}_sesid.txt | awk {'print $1'} | sed "${s}q;d")


                                        echo "----"
                                        echo ${TempTreMethod}, ${METH_ID}, ${ARO}, ${s}, ${SubID}


                                        SubDIR="${ConDir}/${SubID}_${SesID}"

                                        TDIR=${SubDIR}/ED${EDtype}_${METH_ID}_AR${ARO}_MA${MAO}_FWHM${FWHMsize}_${TempTreMethod}${TempTrendP[$tr_cnt]}_cbhat_ICACLEAN0_GSR0.nii.gz
					wTDIR=${SubDIR}/ED${EDtype}_${METH_ID}_AR${ARO}_MA${MAO}_FWHM${FWHMsize}_${TempTreMethod}${TempTrendP[$tr_cnt]}_Wcbhat_ICACLEAN0_GSR0.nii.gz
					wVDIR=${SubDIR}/ED${EDtype}_${METH_ID}_AR${ARO}_MA${MAO}_FWHM${FWHMsize}_${TempTreMethod}${TempTrendP[$tr_cnt]}_wvc_ICACLEAN0_GSR0.nii.gz

					if [ $METH_ID == '3dREMLfit' ]; then
						echo "copy stuff to the right naming format."
						WBETAtmp=${SubDIR}/${EDtype}_sub-${SubID}_ses-${SesID}_Decon_wbeta.nii.gz
						WSTDtmp=${SubDIR}/${EDtype}_sub-${SubID}_ses-${SesID}_Decon_REMLvar_wstd.nii.gz

						${FSLDIR}/bin/fslmaths ${WBETAtmp} ${wTDIR} #just copy
						${FSLDIR}/bin/fslmaths ${WSTDtmp} -sqr ${wVDIR} # AFNI outputs std, we need variance
 					fi

#                                        if [ ! -f $TDIR ] || [ ! -f $wTDIR ] || [ ! -f $wVDIR ]; then
					if  [ ! -f $wTDIR ] || [ ! -f $wVDIR ]; then
                                                echo "$TDIR does not exists."
                                                continue 1
                                        fi

					if [ $COHORT == "ROCKLAND" ]; then
						InDir="/well/nichols/users/scf915/${COHORT}/R_mpp/sub-${SubID}/ses-${SesID}/sub-${SubID}_ses-${SesID}_${EDtype_Rmpp}_bold_mpp"
					elif [ $COHORT == "tHCP" ]; then
						InDir="/well/nichols/users/scf915/${COHORT}/R_mpp/sub-${SubID}/ses-${SesID}/${SubID}_3T_tfMRI_${SesID}_mpp" #100307_3T_tfMRI_MOTOR_LR_mpp
					fi
#					echo $InDir

#					sh func2mni.sh $TDIR $InDir # We don't really need the none-whitened case
					sh func2mni.sh $wTDIR $InDir
					sh func2mni.sh $wVDIR $InDir


                                echo "---"
                        done

			#fslmerge for each scenario
			LEVEL2DIR=${ConDir}/LEVEL2

			mkdir -p ${LEVEL2DIR}


			# We don't really need the non-whitened case! AND even if we need, only a loop around detrending methods should be enough.
			#${FSLDIR}/bin/fslmerge \
			#-t  ${LEVEL2DIR}/ED${EDtype}_${METH_ID}_AR${ARO}_MA${MAO}_FWHM${FWHMsize}_${TempTreMethod}${TempTrendP[$tr_cnt]}_cbhat_4D \
			#${ConDir}/*_*/mni/*_cbhat_ICACLEAN0_GSR0_mni.nii.gz
			#echo "4D data is: ${LEVEL2DIR}/ED${EDtype}_${METH_ID}_AR${ARO}_MA${MAO}_FWHM${FWHMsize}_${TempTreMethod}${TempTrendP[$tr_cnt]}_cbhat_4D"

			${FSLDIR}/bin/fslmerge \
			-t  ${LEVEL2DIR}/ED${EDtype}_${METH_ID}_AR${ARO}_MA${MAO}_FWHM${FWHMsize}_${TempTreMethod}${TempTrendP[$tr_cnt]}_Wcbhat_4D \
			${ConDir}/*_*/mni/*_Wcbhat_ICACLEAN0_GSR0_mni.nii.gz
                        echo "4D data is: ${LEVEL2DIR}/ED${EDtype}_${METH_ID}_AR${ARO}_MA${MAO}_FWHM${FWHMsize}_${TempTreMethod}${TempTrendP[$tr_cnt]}_Wcbhat_4D"

                        ${FSLDIR}/bin/fslmerge \
                        -t  ${LEVEL2DIR}/ED${EDtype}_${METH_ID}_AR${ARO}_MA${MAO}_FWHM${FWHMsize}_${TempTreMethod}${TempTrendP[$tr_cnt]}_wvc_4D \
                        ${ConDir}/*_*/mni/*_wvc_ICACLEAN0_GSR0_mni.nii.gz
                        echo "4D data is: ${LEVEL2DIR}/ED${EDtype}_${METH_ID}_AR${ARO}_MA${MAO}_FWHM${FWHMsize}_${TempTreMethod}${TempTrendP[$tr_cnt]}_wvc_4D"




                done
        done
        tr_cnt=$((tr_cnt + 1))
done

