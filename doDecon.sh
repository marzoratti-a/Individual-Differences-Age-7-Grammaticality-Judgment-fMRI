#!/bin/tcsh

# first cd /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/upload

if ( $#argv > 0 ) then
    set subj = $argv[1]
else
    set subj = s01
endif

3dDeconvolve -input r*_scale.nii                            \
#    -censor motion_${subj}_censor.1D                                         \
    -mask full_mask.nii						     \
    -polort 2                                                                \
    -num_stimts 14                                                           \
    -stim_times 1 stimuli/sem.1D 'BLOCK(2,1)'                          \
    -stim_label 1 sem                                                  \
    -stim_times 2 stimuli/plaus.1D 'BLOCK(2,1)'                        \
    -stim_label 2 plaus                                                \
    -stim_times 3 stimuli/gram.1D 'BLOCK(2,1)'                        \
    -stim_label 3 gram                                                \
    -stim_file 4 regressors/trans_x_sem_all_tmp.txt'[0]' -stim_base 4 -stim_label 4 trans_x_sem   \
    -stim_file 5 regressors/trans_y_sem_all_tmp.txt'[0]' -stim_base 5 -stim_label 5 trans_y_sem   \
    -stim_file 6 regressors/trans_z_sem_all_tmp.txt'[0]' -stim_base 6 -stim_label 6 trans_z_sem   \
    -stim_file 7 regressors/rot_x_sem_all_tmp.txt'[0]' -stim_base 7 -stim_label 7 rot_x_sem     \
    -stim_file 8 regressors/rot_y_sem_all_tmp.txt'[0]' -stim_base 8 -stim_label 8 rot_y_sem     \
    -stim_file 9 regressors/rot_z_sem_all_tmp.txt'[0]' -stim_base 9 -stim_label 9 rot_z_sem     \
    -stim_file 10 regressors/global_signal_sem_all_tmp.txt'[0]' -stim_base 10 -stim_label 10 global_signal_sem     \
    -stim_file 11 regressors/csf_sem_all_tmp.txt'[0]' -stim_base 11 -stim_label 11 csf_sem     \
    -stim_file 12 regressors/white_matter_sem_all_tmp.txt'[0]' -stim_base 12 -stim_label 12 white_matter_sem     \
    -stim_file 13 regressors/trans_x_plaus_all_tmp.txt'[0]' -stim_base 13 -stim_label 13 trans_x_plaus   \
    -stim_file 14 regressors/trans_y_plaus_all_tmp.txt'[0]' -stim_base 14 -stim_label 14 trans_y_plaus    \
    -stim_file 15 regressors/trans_z_plaus_all_tmp.txt'[0]' -stim_base 15 -stim_label 15 trans_z_plaus \
    -stim_file 16 regressors/rot_x_plaus_all_tmp.txt'[0]' -stim_base 16 -stim_label 16 rot_x_plaus  \
    -stim_file 17 regressors/rot_y_plaus_all_tmp.txt'[0]' -stim_base 17 -stim_label 17 rot_y_plaus  \
    -stim_file 18 regressors/rot_z_plaus_all_tmp.txt'[0]' -stim_base 18 -stim_label 18 rot_z_plaus  \
    -stim_file 19 regressors/global_signal_plaus_all_tmp.txt'[0]' -stim_base 19 -stim_label 19 global_signal_plaus     \
    -stim_file 20 regressors/csf_plaus_all_tmp.txt'[0]' -stim_base 20 -stim_label 20 csf_plaus     \
    -stim_file 21 regressors/white_matter_plaus_all_tmp.txt'[0]' -stim_base 21 -stim_label 21 white_matter_plaus     \
    -stim_file 22 regressors/trans_x_gram_all_tmp.txt'[0]' -stim_base 22 -stim_label 22 trans_x_gram   \
    -stim_file 23 regressors/trans_y_gram_all_tmp.txt'[0]' -stim_base 23 -stim_label 23 trans_y_gram    \
    -stim_file 24 regressors/trans_z_gram_all_tmp.txt'[0]' -stim_base 24 -stim_label 24 trans_z_gram \
    -stim_file 25 regressors/rot_x_gram_all_tmp.txt'[0]' -stim_base 25 -stim_label 25 rot_x_gram  \
    -stim_file 26 regressors/rot_y_gram_all_tmp.txt'[0]' -stim_base 26 -stim_label 26 rot_y_gram  \
    -stim_file 27 regressors/rot_z_gram_all_tmp.txt'[0]' -stim_base 27 -stim_label 27 rot_z_gram  \
    -stim_file 28 regressors/global_signal_gram_all_tmp.txt'[0]' -stim_base 28 -stim_label 28 global_signal_gram     \
    -stim_file 29 regressors/csf_gram_all_tmp.txt'[0]' -stim_base 29 -stim_label 29 csf_gram     \
    -stim_file 30 regressors/white_matter_gram_all_tmp.txt'[0]' -stim_base 30 -stim_label 30 white_matter_gram     \
    -jobs 8                                                                  \
    -gltsym 'SYM: sem -plaus'				     \
    -glt_label 1 sem -plaus					     \
    -gltsym 'SYM: plaus -gram'				     \
    -glt_label 2 plaus-gram					     \
    -gltsym 'SYM: sem -gram'				     \
    -glt_label 3 sem-gram					     \
    -fout -tout -x1D X.xmat.1D -xjpeg X.jpg                                  \
    -x1D_uncensored X.nocensor.xmat.1D                                       \
    -fitts fitts.$subj                                                       \
    -errts errts.${subj}                                                     \
    -bucket stats.$subj
