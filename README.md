# fMRI-Scripts-DMPM

1. Copy event files from LabWork folder into func folder & rename regressor files to be {$task}_regressors.tsv
   Copy bold mask files into stimuli folder, rename to ${task}_mask.nii.gz

2. Run script to change time files to AFNI format
   cd /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/upload
   bash make_Timings.sh

3. Run smoothing and scaling script (3dmerge, 3dTstat, 3dcalc)

4. Run script to get confound .txt files
   cd /mnt/c/users/anama/onedrive/documents/labwork/proj-62bddf5ef3194eded6f9293d/bids/derivatives/upload
   bash get_Confounds.sh
   
5. Run script to complete first level analysis/comparisons
   doDecon.sh
