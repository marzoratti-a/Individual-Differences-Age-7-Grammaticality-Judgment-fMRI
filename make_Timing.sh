#!/bin/bash


#Check whether the file subjList.txt exists; if not, create it
if [ ! -f subjList.txt ]; then
	ls | grep ^sub- > subjList.txt
fi

#Loop over all subjects and format timing files into FSL format
for subj in `cat subjList.txt`; do
	cd /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f92s/upload/$subj/ses-7/func
	cat ${subj}_ses-7_task-Sem_* | awk '{if ($3=="S_C") {print $1, $2, 1}}' > sem_ctrl_run1.txt
	cat ${subj}_ses-7_task-Sem_* | awk '{if ($3=="S_H") {print $1, $2, 1}}' > sem_h_run1.txt
	cat ${subj}_ses-7_task-Sem_* | awk '{if ($3=="S_L") {print $1, $2, 1}}' > sem_l_run1.txt
	cat ${subj}_ses-7_task-Sem_* | awk '{if ($3=="S_U") {print $1, $2, 1}}' > sem_u_run1.txt
	
	cat ${subj}_ses-7_task-Plaus_* | awk '{if ($3=="SP_C") {print $1, $2, 1}}' > plaus_ctrl_run1.txt
	cat ${subj}_ses-7_task-Plaus_* | awk '{if ($3=="SP_I") {print $1, $2, 1}}' > plaus_i_run1.txt
	cat ${subj}_ses-7_task-Plaus_* | awk '{if ($3=="SP_S") {print $1, $2, 1}}' > plaus_s_run1.txt
	cat ${subj}_ses-7_task-Plaus_* | awk '{if ($3=="SP_W") {print $1, $2, 1}}' > plaus_w_run1.txt
	
	cat ${subj}_ses-7_task-Gram_* | awk '{if ($3=="G_C") {print $1, $2, 1}}' > gram_ctrl_run1.txt
	cat ${subj}_ses-7_task-Gram_* | awk '{if ($3=="G_F") {print $1, $2, 1}}' > gram_f_run1.txt
	cat ${subj}_ses-7_task-Gram_* | awk '{if ($3=="G_G") {print $1, $2, 1}}' > gram_g_run1.txt
	cat ${subj}_ses-7_task-Gram_* | awk '{if ($3=="G_P") {print $1, $2, 1}}' > gram_p_run1.txt
	
#Now convert to AFNI format
	timing_tool.py -fsl_timing_files sem*.txt -write_timing sem.1D
	timing_tool.py -fsl_timing_files plaus*.txt -write_timing plaus.1D
	timing_tool.py -fsl_timing_files plaus*.txt -write_timing gram.1D

	cd ../..

done
