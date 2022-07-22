#!/bin/tcsh

# first cd /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/upload
for subj in `cat subjList.txt`; do
cd /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/output/$subj/ses-7/func

for task in sem plaus gram; do

if ( $#argv > 0 ) then
    set subj = $argv[1]
else
    set subj = s01
endif

3dDeconvolve -input r*_scale.nii                            \
#    -censor motion_${subj}_censor.1D                                         \
    -mask full_mask.nii						     \
    -polort 1                                                                \
    -num_stimts 30                                                           \
    -fout                                       \
    -rout                                       \
    -stim_times 1 stimuli/sem.1D 'BLOCK(2,1)'                          \
    -stim_label 1 sem                                                  \
    -stim_times 2 stimuli/sem_ctrl.1D 'BLOCK(2,1)'                        \
    -stim_label 2 sem_ctrl                                                \
    -stim_file 4 regressors/trans_x_${task}_all_tmp.txt'[0]' -stim_base 4 -stim_label 4 trans_x_${task}   \
    -stim_file 5 regressors/trans_y_${task}_all_tmp.txt'[0]' -stim_base 5 -stim_label 5 trans_y_${task}   \
    -stim_file 6 regressors/trans_z_${task}_all_tmp.txt'[0]' -stim_base 6 -stim_label 6 trans_z_${task}   \
    -stim_file 7 regressors/rot_x_${task}_all_tmp.txt'[0]' -stim_base 7 -stim_label 7 rot_x_${task}     \
    -stim_file 8 regressors/rot_y_${task}_all_tmp.txt'[0]' -stim_base 8 -stim_label 8 rot_y_${task}    \
    -stim_file 9 regressors/rot_z_${task}_all_tmp.txt'[0]' -stim_base 9 -stim_label 9 rot_z_${task}    \
    -stim_file 10 regressors/global_signal_${task}_all_tmp.txt'[0]' -stim_base 10 -stim_label 10 global_signal_${task}     \
    -stim_file 11 regressors/csf_${task}_all_tmp.txt'[0]' -stim_base 11 -stim_label 11 csf_sem     \
    -stim_file 12 regressors/white_matter_${task}_all_tmp.txt'[0]' -stim_base 12 -stim_label 12 white_matter_${task}     \
    -jobs 8                                                                  \
    -gltsym 'SYM: ${task} -sem_ctrl'				     \
    -glt_label 1 ${task} -ctrl					     \
    -gltsym 'SYM: ${task}_ctrl -sem'				     \
    -glt_label 2 ctrl -sem					     \
  
    -fout -tout -x1D X.xmat.1D -xjpeg X.jpg                                  \
    -x1D_uncensored X.nocensor.xmat.1D                                       \
    -fitts fitts.$subj                                                       \
    -errts errts.${subj}                                                     \
    -bucket stats.$subj

    done
done
