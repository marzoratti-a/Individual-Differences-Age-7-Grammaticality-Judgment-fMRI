#!/bin/bash

# first cd /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d

#Check whether the file subjList.txt exists; if not, create it
if [ ! -f subjList.txt ]; then
	ls | grep ^sub- > subjList.txt
fi

for subj in `cat subjList.txt`; do
  cd /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/$subj
  subjid=${subj%.ses-7*}
  
   if [ ! -d bold ]; then
      mkdir bold
   fi
   if [ ! -f bold.txt ]; then
      ls | grep ^dt-neuro-func-task > bold.txt
   fi

    for bold_folder in `cat bold.txt`; do
      cd ${bold_folder}
      task=$(
          awk '
              match($0,/Sem|Plaus|Gram/) {
                  print tolower(substr($0,RSTART,RLENGTH))
                  exit
              }
          ' "_info.json"
      )
      echo $task $subj
      mv ".nii.gz" "$subj_$task_bold.nii.gz" 
      cp "$subj_$task_bold.nii.gz" ../bold
      cd ..
     done

     cd bold
     
     cp -r /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/$subj/bold /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/upload/$subjid/ses-7/func
   
    done
