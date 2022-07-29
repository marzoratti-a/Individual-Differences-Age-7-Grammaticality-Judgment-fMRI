#!/bin/tcsh

# first go to subject folder 
# cd /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/upload/sub-

    if ( $#argv > 0 ) then
        set subj = $argv[1]
    else
        set subj = s01
    endif
    
 3dDeconvolve -input rsem_scale.nii                            \
        -mask stimuli/sem_full_mask.nii.gz					     \
        -polort 1                                                                \
        -xout -progress                                                     \
        -num_stimts 11                                                           \
        -stim_times 1 stimuli/sem.1D 'BLOCK(2,1)'                          \
        -stim_label 1 sem                                                  \
        -stim_times 2 stimuli/sem_ctrl.1D 'BLOCK(2,1)'                          \
        -stim_label 2 sem_ctrl                                               \
        -stim_file 3 regressors/trans_x_sem_noHead_tmp.txt'[0]' -stim_base 3 -stim_label 3 trans_x_sem   \
        -stim_file 4 regressors/trans_y_sem_noHead_tmp.txt'[0]' -stim_base 4 -stim_label 4 trans_y_sem   \
        -stim_file 5 regressors/trans_z_sem_noHead_tmp.txt'[0]' -stim_base 5 -stim_label 5 trans_z_sem   \
        -stim_file 6 regressors/rot_x_sem_noHead_tmp.txt'[0]' -stim_base 6 -stim_label 6 rot_x_sem     \
        -stim_file 7 regressors/rot_y_sem_noHead_tmp.txt'[0]' -stim_base 7 -stim_label 7 rot_y_sem     \
        -stim_file 8 regressors/rot_z_sem_noHead_tmp.txt'[0]' -stim_base 8 -stim_label 8 rot_z_sem     \
        -stim_file 9 regressors/global_signal_sem_noHead_tmp.txt'[0]' -stim_base 9 -stim_label 9 global_signal_sem     \
        -stim_file 10 regressors/csf_sem_noHead_tmp.txt'[0]' -stim_base 10 -stim_label 10 csf_sem     \
        -stim_file 11 regressors/white_matter_sem_noHead_tmp.txt'[0]' -stim_base 11 -stim_label 11 white_matter_sem     \
        -jobs 8  \
        # tells AFNI to run GLT based on label names
        -gltsym 'SYM: sem -sem_ctrl'				     \
        -glt_label 1 sem -ctrl					     \
        -gltsym 'SYM: sem_ctrl -sem'				     \
        -glt_label 2 ctrl -sem					     \

        -fout -tout -rout -x1D X.xmat.1D -xjpeg X.jpg                                  \
        -x1D_uncensored X.nocensor.xmat.1D                                       \
        -fitts fitts.$subj.sem                                                       \
        -errts errts.${subj}.sem                                                     \
        -bucket stats.$subj.sem                                                   \

 3dDeconvolve -input rplaus_scale.nii                            \
        -mask stimuli/plaus_full_mask.nii.gz					     \
        -polort 1                                                                \
        -xout -progress                                                     \
        -num_stimts 11                                                           \
        -stim_times 1 stimuli/plaus.1D 'BLOCK(2,1)'                          \
        -stim_label 1 plaus                                                  \
        -stim_times 2 stimuli/plaus_ctrl.1D 'BLOCK(2,1)'                          \
        -stim_label 2 plaus_ctrl                                               \
        -stim_file 3 regressors/trans_x_plaus_noHead_tmp.txt'[0]' -stim_base 3 -stim_label 3 trans_x_plaus   \
        -stim_file 4 regressors/trans_y_plaus_noHead_tmp.txt'[0]' -stim_base 4 -stim_label 4 trans_y_plaus   \
        -stim_file 5 regressors/trans_z_plaus_noHead_tmp.txt'[0]' -stim_base 5 -stim_label 5 trans_z_plaus   \
        -stim_file 6 regressors/rot_x_plaus_noHead_tmp.txt'[0]' -stim_base 6 -stim_label 6 rot_x_plaus     \
        -stim_file 7 regressors/rot_y_plaus_noHead_tmp.txt'[0]' -stim_base 7 -stim_label 7 rot_y_plaus     \
        -stim_file 8 regressors/rot_z_plaus_noHead_tmp.txt'[0]' -stim_base 8 -stim_label 8 rot_z_plaus     \
        -stim_file 9 regressors/global_signal_plaus_noHead_tmp.txt'[0]' -stim_base 9 -stim_label 9 global_signal_plaus     \
        -stim_file 10 regressors/csf_plaus_noHead_tmp.txt'[0]' -stim_base 10 -stim_label 10 csf_plaus     \
        -stim_file 11 regressors/white_matter_plaus_noHead_tmp.txt'[0]' -stim_base 11 -stim_label 11 white_matter_plaus     \
        -jobs 8  \
        # tells AFNI to run GLT based on label names
        -gltsym 'SYM: plaus -plaus_ctrl'				     \
        -glt_label 1 plaus -ctrl					     \
        -gltsym 'SYM: plaus_ctrl -plaus'				     \
        -glt_label 2 ctrl -plaus					     \

        -fout -tout -rout -x1D X.xmat.1D -xjpeg X.jpg                                  \
        -x1D_uncensored X.nocensor.xmat.1D                                       \
        -fitts fitts.$subj.plaus                                                      \
        -errts errts.${subj}.plaus                                                    \
        -bucket stats.$subj.plaus                        \
        
  3dDeconvolve -input rgram_scale.nii                            \
        -mask stimuli/gram_full_mask.nii.gz					     \
        -polort 1                                                                \
        -xout -progress                                                     \
        -num_stimts 11                                                           \
        -stim_times 1 stimuli/gram.1D 'BLOCK(2,1)'                          \
        -stim_label 1 gram                                                  \
        -stim_times 2 stimuli/gram_ctrl.1D 'BLOCK(2,1)'                          \
        -stim_label 2 gram_ctrl                                               \
        -stim_file 3 regressors/trans_x_gram_noHead_tmp.txt'[0]' -stim_base 3 -stim_label 3 trans_x_gram   \
        -stim_file 4 regressors/trans_y_gram_noHead_tmp.txt'[0]' -stim_base 4 -stim_label 4 trans_y_gram   \
        -stim_file 5 regressors/trans_z_gram_noHead_tmp.txt'[0]' -stim_base 5 -stim_label 5 trans_z_gram   \
        -stim_file 6 regressors/rot_x_gram_noHead_tmp.txt'[0]' -stim_base 6 -stim_label 6 rot_x_gram     \
        -stim_file 7 regressors/rot_y_gram_noHead_tmp.txt'[0]' -stim_base 7 -stim_label 7 rot_y_gram     \
        -stim_file 8 regressors/rot_z_gram_noHead_tmp.txt'[0]' -stim_base 8 -stim_label 8 rot_z_gram     \
        -stim_file 9 regressors/global_signal_gram_noHead_tmp.txt'[0]' -stim_base 9 -stim_label 9 global_signal_gram     \
        -stim_file 10 regressors/csf_gram_noHead_tmp.txt'[0]' -stim_base 10 -stim_label 10 csf_gram     \
        -stim_file 11 regressors/white_matter_gram_noHead_tmp.txt'[0]' -stim_base 11 -stim_label 11 white_matter_gram     \
        -jobs 8  \
        # tells AFNI to run GLT based on label names
        -gltsym 'SYM: gram -gram_ctrl'				     \
        -glt_label 1 gram -ctrl					     \
        -gltsym 'SYM: gram_ctrl -gram'				     \
        -glt_label 2 ctrl -gram					     \

        -fout -tout -rout -x1D X.xmat.1D -xjpeg X.jpg                                  \
        -x1D_uncensored X.nocensor.xmat.1D                                       \
        -fitts fitts.$subj.gram                                                       \
        -errts errts.${subj}.gram                                                     \
        -bucket stats.$subj.gram                        \
