#!/bin/bash

# first cd /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/upload
# Check whether the file subjList.txt exists; if not, create it
if [ ! -f subjList.txt ]; then
	ls | grep ^sub- > subjList.txt
fi

for subj in `cat subjList.txt`; do
  cd /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/brainlife.app-fmriprep/$subj/ses-7/notspecific

find . -type f -name "*_regressors.json" -print0 | while IFS= read -r -d '' f; do
    if str=$(grep -wE "Sem|Plaus|Gram" "$f"); then              # search the json file for the strings
        str=$(head -n 1 <<< "$str" | tr [:upper:] [:lower:])    # pick the 1st match and lower the case
        base=${f%.json}                                         # remove the extention
        echo mv -- "${base}.tsv" "${base}_${str}.tsv"           # rename the file
    fi
done
  
  
