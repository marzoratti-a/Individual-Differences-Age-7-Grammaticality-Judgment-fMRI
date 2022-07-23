#!/bin/tcsh

# first cd /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/upload
for subj in `cat subjList.txt`; do
cd /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/upload/$subj/ses-7/func

if ( $#argv > 0 ) then
    set subj = $argv[1]
else
    set subj = s01
endif

for task in sem plaus gram; do

3dDeconvolve -input *_tag-${task}_*_bold.nii.gz                            \
#    -censor motion_${subj}_censor.1D                                         \
    -mask ${task}_mask.nii.gz					     \
    -polort 1                                                                \
    -xout -progress                                                     \
    -num_stimts 11                                                           \
    -stim_times 1 stimuli/${task}.1D 'BLOCK(2,1)'                          \
    -stim_label 1 ${task}                                                  \
    -stim_times 2 stimuli/${task}_ctrl.1D 'BLOCK(2,1)'                        \
    -stim_label 2 ${task}_ctrl                                                \
    -stim_file 3 regressors/trans_x_${task}_noHead_tmp.txt'[0]' -stim_base 3 -stim_label 3 trans_x_${task}   \
    -stim_file 4 regressors/trans_y_${task}_noHead_tmp.txt'[0]' -stim_base 4 -stim_label 4 trans_y_${task}   \
    -stim_file 5 regressors/trans_z_${task}_noHead_tmp.txt'[0]' -stim_base 5 -stim_label 5 trans_z_${task}   \
    -stim_file 6 regressors/rot_x_${task}_noHead_tmp.txt'[0]' -stim_base 6 -stim_label 6 rot_x_${task}     \
    -stim_file 7 regressors/rot_y_${task}_noHead_tmp.txt'[0]' -stim_base 7 -stim_label 7 rot_y_${task}    \
    -stim_file 8 regressors/rot_z_${task}_noHead_tmp.txt'[0]' -stim_base 8 -stim_label 8 rot_z_${task}    \
    -stim_file 9 regressors/global_signal_${task}_noHead_tmp.txt'[0]' -stim_base 9 -stim_label 9 global_signal_${task}     \
    -stim_file 10 regressors/csf_${task}_noHead_tmp.txt'[0]' -stim_base 10 -stim_label 10 csf_sem     \
    -stim_file 11 regressors/white_matter_${task}_noHead_tmp.txt'[0]' -stim_base 11 -stim_label 11 white_matter_${task}     \
    -jobs 8  \
    # tells AFNI to run GLT based on label names
    -gltsym 'SYM: ${task} -${task}_ctrl'				     \
    -glt_label 1 ${task} -ctrl					     \
    -gltsym 'SYM: ${task}_ctrl -${task}'				     \
    -glt_label 2 ctrl -${task}					     \
  
    -fout -tout -rout -x1D X.xmat.1D -xjpeg X.jpg                                  \
    -x1D_uncensored X.nocensor.xmat.1D                                       \
    -fitts fitts.$subj                                                       \
    -errts errts.${subj}                                                     \
    -bucket stats.$subj

    done
done
