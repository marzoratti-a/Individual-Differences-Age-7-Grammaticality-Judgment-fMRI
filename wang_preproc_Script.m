% Script created by Andrew Jahn, University of Michigan, 10 January 2020
% This script is designed to analyze the Flanker dataset (https://openneuro.org/datasets/ds000102/versions/00001)
% It is used during the Scripting chapter of the SPM walkthrough (https://andysbrainbook.readthedocs.io/en/latest/SPM/SPM_Short_Course/SPM_06_Scripting.html)
% -------------------------------------------------------------------------
% SET UP MATLABBATCH
% 1) Raw images should be in func folder within subject directory,
% structural T1w in anat folder
% 2) Realign (Estimate and Reslice)--> rsub...bold.nii (sessions + mean)
% 3) Slice Time Correction--> arsub...bold.nii
% 4) Coregister (Estimate and Reslice)--> rsub...T1w.nii
% 5) Segment the anatomical data for white, grey, csf, soft tissue, bone,
% other--> c#sub...T1w.nii, y_sub...T1w.nii (deformation)
% 6) Normalize slice-time corrected func files to c1, c2, and c3 to the T1
% template (y_sub...T1w.nii)--> warsub.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            nbbbbbbbbbbbbbbbbbbbbbbbbbbbbb gmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmd                                                                           sg..bold.nii
% 7) Smoothing--> swarsub...bold.nii
% 8) Splits Session 1 and 2 files for first level analysis
% 9) Model Specification (defines conditions and timing)
% 10) Model Estimation
% 11) Contrast Manager (defines condition-control contrasts + condition)
% 12) Deletes all except the final swarsub...bold.nii files from func folder
% -------------------------------------------------------------------------
clc;
clear all;
rootdir="D:\agespec_data";
%Set as directory with fMRI subject folders
d=dir(rootdir); 
cd(rootdir);
%spm fmri
% addpath('C:\Program Files\spm12')

for i=[64]
    % 5192 no Sem ses 2
    subject = num2str(d(i).name,'%02d'); % Zero-pads each number so that the subject ID is 2 characters long

%     % List of open inputs
% 
%     %% GRAM TASK Processing
%     nrun = 1; % enter the number of runs here
%     jobfile = {'D:\MATLAB\wang_preprocG_job.m'};
%     spm('defaults', 'FMRI');
%     spm_jobman('initcfg');
%     run(jobfile); %creates matlabbatch cell structure
%     spm_jobman('run', matlabbatch);
%     
%     clear matlabbatch jobfile jobs 
%     
%     %% SEM TASK Processing
%         
%     % List of open inputs
%     nrun = 1; % enter the number of runs xhere
%     jobfile = {'D:\MATLAB\wang_preprocS_job.m'};
%     spm('defaults', 'FMRI');
%     spm_jobman('initcfg');
%     run(jobfile); %creates matlabbatch cell structure
%     spm_jobman('run', matlabbatch);
%     
%     clear matlabbatch jobfile jobs 
    
        
    %% PLAUS TASK Processing
        
    % List of open inputs
    nrun = 1; % enter the number of runs here
    jobfile = {'D:\MATLAB\wang_preprocP_job.m'};
    spm('defaults', 'FMRI');
    spm_jobman('initcfg');
    run(jobfile); %creates matlabbatch cell structure
    spm_jobman('run', matlabbatch);
    
    clear matlabbatch jobfile jobs 
    
      
 
end

disp('All files done, rejoice and drink')