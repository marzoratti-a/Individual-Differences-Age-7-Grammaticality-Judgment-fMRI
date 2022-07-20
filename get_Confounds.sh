#!/bin/bash
# first cd /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/upload
# Check whether the file subjList.txt exists; if not, create it
if [ ! -f subjList.txt ]; then
	ls | grep ^sub- > subjList.txt
fi

for subj in `cat subjList.txt`; do
  cd /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/brainlife.app-fmriprep/$subj/ses-7/notspecific
  
  # Make txt files with and without header for each regressor of interest
  for reg in global_signal csf white_matter trans_x trans_y trans_z rot_x rot_y rot_z; do
	  for task in sem plaus gram; do
	    awk -v col=$reg 'NR==1{for(i=1;i<=NF;i++){if($i==col){c=i;break}} print $c} NR>1{print $c}' *_${task}_regressors.tsv > ${reg}_${task}_all_tmp.txt;
	    sed '1d' ${reg}_${task}_all_tmp.txt > ${reg}_${task}_noHead_tmp.txt
	    sed '1!d' ${reg}_${task}_all_tmp.txt > ${reg}_${task}_Head_tmp.txt
	  done
  done
  
   # Move newly made regressor/confound files to subject folder, rename to regressors
     cp -avr /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/brainlife.app-fmriprep/$subj/ses-7/notspecific /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/upload/$subj
     cd /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/upload/$subj
     mv notspecific regressors
     
     cd ..
  
done
