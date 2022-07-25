#!/bin/bash

# first cd /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/upload
# Check whether the file subjList.txt exists; if not, create it
if [ ! -f subjList.txt ]; then
	ls | grep ^sub- > subjList.txt
fi

for subj in `cat subjList.txt`; do
  cd /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/brainlife.app-fmriprep/$subj/ses-7/notspecific
  
  # https://stackoverflow.com/questions/47056022/how-to-show-filename-when-awk-found-matches
  grep -ro 'Gram' --include=*.json
  $ awk '/Sem/' *_regressors.tsv > *_sem_regressors.tsv
  
for match in "$(grep -ro 'Sem')";do
    echo mv "${match%:*.tsv}" "${match#*:}"
done

# within folder, finds files with json extension and saves name of file as "filename"
# loops through and if file contains 'Sem', adds "_sem" to corresponding file with filename .tsv 

find -type f -name "*_regressors.json" -print0 | while IFS= read -r -d '' filename
do
    if [[grep -q 'Sem' "$filename"]]; then
        sem_name="${filename%.*}" 
	mv ${sem_name}.tsv ${sem_name}_sem.tsv
    fi 
    
    if [[grep -q 'Plaus' "$filename"]]; then
    plaus_name="${filename%.*}"
    mv ${plaus_name}.tsv ${plaus_name}_plaus.tsv
    fi
    
    if [[grep -q 'Gram' "$filename"]]; then
        gram_name="${filename%.*}"
	mv ${gram_name}.tsv ${gram_name}_gram.tsv
    fi
done
  
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
  
done
