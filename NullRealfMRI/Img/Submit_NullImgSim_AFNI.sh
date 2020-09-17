#bin/bash
set -e

AFNI_SPM_singularity_image=/apps/singularity/afni-r-python3-2020-03-26-v1.sif
AFNI_bin=/opt/afni-latest

BOLDIMG=${path_data}/${subject}_${task}_bold.nii
stim_times=
HRF=
smoothing=
zero_aux=
TR=

singularity exec --cleanenv -B $home_dir $AFNI_SPM_singularity_image \
$AFNI_bin/afni_proc.py 							 \
-subj_id $subject                                                        \
-script proc.$subject -scr_overwrite                                     \
-regress_opts_3dD -force_TR $TR                                          \
-tcat_remove_first_trs 0                                                 \
-dsets ${BOLDIMG}			                                 \
-volreg_align_to third                                                   \
-regress_polort 2                                                        \
-regress_bandpass 0.01 10                                                \
-blur_size ${smoothing}.${zero_aux}                                      \
-regress_stim_types AM1                                                  \
-regress_stim_times $stim_times                                          \
-regress_stim_labels activation_stimulus                                 \
-regress_basis ${HRF}                                                    \
-regress_make_ideal_sum sum_ideal.1D                                     \
-regress_run_clustsim no                                                 \
-regress_est_blur_epits                                                  \
-regress_est_blur_errts                                                  \
-regress_reml_exec                                                       \
-regress_opts_reml -Rwherr whitened_errts.${subject}_REML

