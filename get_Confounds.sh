#!/bin/bash

cd /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/upload

#Check whether the file subjList.txt exists; if not, create it
if [ ! -f subjList.txt ]; then
	ls | grep ^sub- > subjList.txt
fi

for subj in `cat subjList.txt`; do
	cd $subj
  # Make backup versions of regressor/confound files
  cp /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/brainlife.app-fmriprep/$subj/ses-7/notspecific/*sem_regressors.tsv /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/upload/$subj/stimuli/sem_regressors.tsv 
  cp /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/brainlife.app-fmriprep/$subj/ses-7/notspecific/*plaus_regressors.tsv /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/upload/$subj/stimuli/plaus_regressors.tsv 
  cp /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/brainlife.app-fmriprep/$subj/ses-7/notspecific/*gram_regressors.tsv /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/upload/$subj/stimuli/gram_regressors.tsv
  
 for reg in global_signal csf white_matter trans_x trans_y trans_z rot_x rot_y rot_z; do
  for run in 1 2; do
    awk -v col=$reg 'NR==1{for(i=1;i<=NF;i++){if($i==col){c=i;break}} print $c} NR>1{print $c}' sub-08_task-flanker_run-${run}_desc-confounds_regressors.tsv > ${reg}_run${run}_all_tmp.txt;
    sed '1d' ${reg}_run${run}_all_tmp.txt > ${reg}_run${run}_noHead_tmp.txt
    sed '1!d' ${reg}_run${run}_all_tmp.txt > ${reg}_run${run}_Head_tmp.txt
  done
done

done
