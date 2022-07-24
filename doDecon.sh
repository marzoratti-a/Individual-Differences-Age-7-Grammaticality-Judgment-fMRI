#!/bin/tcsh

# first go to subject folder 
# cd /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/upload/sub-

    if ( $#argv > 0 ) then
        set subj = $argv[1]
    else
        set subj = s01
    endif
    
    
foreach task (sem plaus gram)
 3dDeconvolve -input r${task}_scale.nii                            \
        -mask stimuli/${task}_mask.nii.gz					     \
        -polort 1                                                                \
        -xout -progress                                                     \
        -num_stimts 15                                                           \
        -stim_times 1 stimuli/${task}.1D 'BLOCK(2,1)'                          \
        -stim_label 1 ${task}                                               \
        -stim_times 2 stimuli/${task}_ctrl.1D 'BLOCK(2,1)'                          \
        -stim_label 2 ${task}_ctrl       \
        -stim_file 7 regressors/trans_x_${task}_noHead_tmp.txt'[0]' -stim_base 4 -stim_label 4 trans_x_${task}   \
        -stim_file 8 regressors/trans_y_${task}_noHead_tmp.txt'[0]' -stim_base 5 -stim_label 5 trans_y_${task}   \
        -stim_file 9 regressors/trans_z_${task}_noHead_tmp.txt'[0]' -stim_base 6 -stim_label 6 trans_z_${task}   \
        -stim_file 10 regressors/rot_x_${task}_noHead_tmp.txt'[0]' -stim_base 7 -stim_label 7 rot_x_${task}     \
        -stim_file 11 regressors/rot_y_${task}_noHead_tmp.txt'[0]' -stim_base 8 -stim_label 8 rot_y_${task}     \
        -stim_file 12 regressors/rot_z_${task}_noHead_tmp.txt'[0]' -stim_base 9 -stim_label 9 rot_z_${task}    \
        -stim_file 13 regressors/global_signal_${task}_noHead_tmp.txt'[0]' -stim_base 10 -stim_label 10 global_signal_${task}     \
        -stim_file 14 regressors/csf_${task}_noHead_tmp.txt'[0]' -stim_base 11 -stim_label 11 csf_${task}     \
        -stim_file 15 regressors/white_matter_${task}_noHead_tmp.txt'[0]' -stim_base 12 -stim_label 12 white_matter_${task}     \
        -jobs 8  \
        # tells AFNI to run GLT based on label names
        -gltsym 'SYM: ${task} -${task}_ctrl'				     \
        -glt_label 1 ${task} -ctrl					     \
        -gltsym 'SYM: ${task}_ctrl -${task}'				     \
        -glt_label 2 ctrl -${task}					     \
        
        -fout -tout -rout -x1D X.xmat.1D -xjpeg X.jpg                                  \
        -x1D_uncensored X.nocensor.xmat.1D                                       \
        -fitts fitts.${subj}.${task}                                                       \
        -errts errts.${subj}.${task}                                                     \
        -bucket stats.${subj}.${task}                                                     \
end
