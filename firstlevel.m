subjects = [01 02]; % Replace with a list of all of the subjects you wish to analyze

user = getenv('USER'); % Will return the username for OSX operating systems

for subject=subjects

subject = num2str(subject, '%02d');

if exist(['/Users/' user '/Desktop/Flanker/sub-' subject '/func/sub-' subject '_task-flanker_run-1_bold.nii']) == 0
    display('Run 1 has not been unzipped; unzipping now')
    gunzip(['/Users/' user '/Desktop/Flanker/sub-' subject '/func/sub-' subject '_task-flanker_run-1_bold.nii.gz'])
else
    display('Run 1 is already unzipped; doing nothing')
end

if exist(['/Users/' user '/Desktop/Flanker/sub-' subject '/func/sub-' subject '_task-flanker_run-2_bold.nii']) == 0
    display('Run 2 has not been unzipped; unzipping now')
    gunzip(['/Users/' user '/Desktop/Flanker/sub-' subject '/func/sub-' subject '_task-flanker_run-2_bold.nii.gz'])
else
    display('Run 2 is already unzipped; doing nothing')
end

if exist(['/Users/' user '/Desktop/Flanker/sub-' subject '/anat/sub-' subject '_T1w.nii']) == 0
    display('Anatomical image has not been unzipped; unzipping now')
    gunzip(['/Users/' user '/Desktop/Flanker/sub-' subject '/anat/sub-' subject '_T1w.nii.gz'])
else
    display('Anatomical image is already unzipped; doing nothing')
end