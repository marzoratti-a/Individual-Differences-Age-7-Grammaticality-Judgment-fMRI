#!/bin/bash

# first cd /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/upload

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
  cp -avr /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/brainlife.app-fmriprep/$subj/ses-7/notspecific /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/upload/$subj/ses-7/func
  cd /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/upload/$subj/ses-7/func
  mv notspecific regressors
 
  # Loop through tasks for smoothing
  for task in sem plaus gram; do
    3dmerge -1blur_fwhm 6.0 -doall -prefix r${task}_blur.nii \
            ${subj}_*_tag-${task}_*_bold.nii.gz
  done
  
  
  # Loop through tasks for scaling
  for task in sem plaus gram; do
  3dTstat -prefix rm.mean_r${task}.nii r${task}_blur.nii
  3dcalc -a r${task}_blur.nii -b rm.mean_r${task}.nii \
         -c ${subj}_*_tag-${task}_*_bold.nii.gz                          \
         -expr 'c * min(200, a/b*100)*step(a)*step(b)'       \
         -prefix r${task}_scale.nii
  done
  
  rm rm*
  cd ../../..
done
