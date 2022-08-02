#!/bin/bash

cd /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/brainlife.app-fmriprep/sub-5036/ses-7/notspecific

find . -type f -name "*_regressors.json" -print0 |
while IFS= read -r -d '' filename
do
    prefix=${filename%.*}
    suffix=$(
        awk '
            match($0,/Sem|Plaus|Gram/) {
                print tolower(substr($0,RSTART,RLENGTH))
                exit
            }
        ' "$filename"
    )
    mv "$prefix.tsv" "${prefix}_$suffix.tsv" 
done
  
  
