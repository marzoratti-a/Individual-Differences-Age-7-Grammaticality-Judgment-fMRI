#!/bin/bash

# Check whether the file subjList.txt exists; if not, create it
if [ ! -f subjList.txt ]; then
	ls | grep ^sub- > subjList.txt
fi

# Check whether regressors folder exists; if not, create it
if [ ! -d "/mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/upload/$subj/regressors" ]; then
	mkdir regressors
fi

for subj in `cat subjList.txt`; do
  cd /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/brainlife.app-fmriprep/$subj/ses-7/notspecific/
  
  # Make backup versions of regressor/confound files
  cp *sem_regressors.tsv /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/upload/$subj/regressors/$subj_ses-7_sem_regressors.tsv 
  cp *plaus_regressors.tsv /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/upload/$subj/regressors/$subj_ses-7_plaus_regressors.tsv 
  cp *gram_regressors.tsv /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/upload/$subj/regressors/$subj_ses-7_gram_regressors.tsv
  
  cd /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/upload/$subj/regressors
  
  # Make txt files with and without header for each regressor of interest
  for reg in global_signal csf white_matter trans_x trans_y trans_z rot_x rot_y rot_z; do
	  for task in sem plaus gram; do
	    awk -v col=$reg 'NR==1{for(i=1;i<=NF;i++){if($i==col){c=i;break}} print $c} NR>1{print $c}' ${subj}_ses-7_${task}_regressors.tsv > ${reg}_${task}_all_tmp.txt;
	    sed '1d' ${reg}_${task}_all_tmp.txt > ${reg}_${task}_noHead_tmp.txt
	    sed '1!d' ${reg}_${task}_all_tmp.txt > ${reg}_${task}_Head_tmp.txt
	  done
  done
  
  cd ../..
  
done
