 subject = num2str(d(i).name,'%02d'); % Zero-pads each number so that the subject ID is 2 characters long
   %%%%%%%%%%
   matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.name = [subject,'_Gramrun1run2'];
   matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {
                                                                        {['D:\ds003604\',subject,'\ses-7\func\',subject,'_ses-7_task-Gram_acq_run-01_bold.nii']}
                                                                        {['D:\ds003604\',subject,'\ses-7\func\',subject,'_ses-7_task-Gram_acq_run-02_bold.nii']}
                                                                        }';
   % MODULE 2: Realign (Estimate & Reslice)
   matlabbatch{2}.spm.spatial.realign.estwrite.data{1}(1) = cfg_dep(['Named File Selector: ',subject,'_Gramrun1run2(1) - Files'], substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
   matlabbatch{2}.spm.spatial.realign.estwrite.data{2}(1) = cfg_dep(['Named File Selector: ',subject,'_Gramrun1run2(2) - Files'], substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{2}));
   % Specifies run 1 and 2 files selected in named file selector module ^^^
   matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
   matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.sep = 4;
   matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.fwhm = 6;
   % Changed to match Qin et al., 2012 & Szaflarski et al., 2006 ^^^
   matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
   % Registers to mean
   matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.interp = 2;
   % 2nd degree B-spline
   matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
   matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.weight = '';
   % No wrap and no weighting
   matlabbatch{2}.spm.spatial.realign.estwrite.roptions.which = [2 1];
   % Writes out all images + mean image resliced
   matlabbatch{2}.spm.spatial.realign.estwrite.roptions.interp = 4;
   % Interpolates 4th degree B-spline
   matlabbatch{2}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
   matlabbatch{2}.spm.spatial.realign.estwrite.roptions.mask = 1;
   % No wrap, but DOES mask images
   matlabbatch{2}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
 % MODULE 3: Slice Timing
   matlabbatch{3}.spm.temporal.st.scans{1}(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rfiles'));
   matlabbatch{3}.spm.temporal.st.scans{2}(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 2)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{2}, '.','rfiles'));
   % Specifies run 1 and 2 files resliced in previous module' rsub...bold.nii
   matlabbatch{3}.spm.temporal.st.nslices = 56;
   matlabbatch{3}.spm.temporal.st.tr = 1.25;
   matlabbatch{3}.spm.temporal.st.ta = 1.25-(1.25\56); %TR-(TR\# slices)
   matlabbatch{3}.spm.temporal.st.so = [1:2:56 2:2:56]; %interleaved from bottom slice (1) up
   matlabbatch{3}.spm.temporal.st.refslice = 1; %ref slice is first
   matlabbatch{3}.spm.temporal.st.prefix = 'a';
 % MODULE 4: Coregister (Estimate & Reslice)
   matlabbatch{4}.spm.spatial.coreg.estwrite.ref(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
   % Selects mean image produced in realignment\reslice module (Module 2)
   matlabbatch{4}.spm.spatial.coreg.estwrite.source = {['D:\ds003604\' subject '\ses-7\anat\',subject,'_ses-7_acq_T1w.nii']};
   % Selects T1w image for subject from anat directory
   matlabbatch{4}.spm.spatial.coreg.estwrite.other = {''};
   matlabbatch{4}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
   % Objective function= normalized mutual info (default)
   matlabbatch{4}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
   % Avg distance between sampled points in mm; increasingly fine from 4-2mm
   matlabbatch{4}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
   % Tolerance for each parameter
   matlabbatch{4}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
   % Guassian smoothing applied to 256x256mm joint histogram
   matlabbatch{4}.spm.spatial.coreg.estwrite.roptions.interp = 4;
   % T1w image interpolated onto func via 4th degree B-Spline
   matlabbatch{4}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
   matlabbatch{4}.spm.spatial.coreg.estwrite.roptions.mask = 0;
   % No wrapping or masking
   matlabbatch{4}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
   % MODULE 5: Segment
   matlabbatch{5}.spm.spatial.preproc.channel.vols = {['D:\ds003604\',subject,'\ses-7\anat\',subject,'_ses-7_acq_T1w.nii,1']};
   % Selects first T1w image
   matlabbatch{5}.spm.spatial.preproc.channel.biasreg = 0.001;
   % Light bias regularization
   matlabbatch{5}.spm.spatial.preproc.channel.biasfwhm = 60;
   % 60 mm cutoff for Full Width Half Max
   matlabbatch{5}.spm.spatial.preproc.channel.write = [0 1];
   % Set to "Save Bias Corrected" files
   % Below code specifies 6 tissue types (gray matter, white matter, CSF, soft tissue, skull, & all other  (air, tumors, etc)
   % Gray Matter
   matlabbatch{5}.spm.spatial.preproc.tissue(1).tpm = {['C:\Program Files\spm12\tpm\TPM.nii,1']};
   matlabbatch{5}.spm.spatial.preproc.tissue(1).ngaus = 1;
   matlabbatch{5}.spm.spatial.preproc.tissue(1).native = [1 0];
   matlabbatch{5}.spm.spatial.preproc.tissue(1).warped = [0 0];
   % White Matter
   matlabbatch{5}.spm.spatial.preproc.tissue(2).tpm = {['C:\Program Files\spm12\tpm\TPM.nii,2']};
   matlabbatch{5}.spm.spatial.preproc.tissue(2).ngaus = 1;
   matlabbatch{5}.spm.spatial.preproc.tissue(2).native = [1 0];
   matlabbatch{5}.spm.spatial.preproc.tissue(2).warped = [0 0];
   % CSF NIHPD
   matlabbatch{5}.spm.spatial.preproc.tissue(3).tpm = {['C:\Program Files\spm12\tpm\TPM.nii,3']};
   matlabbatch{5}.spm.spatial.preproc.tissue(3).ngaus = 2;
   matlabbatch{5}.spm.spatial.preproc.tissue(3).native = [1 0];
   matlabbatch{5}.spm.spatial.preproc.tissue(3).warped = [0 0];
   % Bone (from SPM)
   matlabbatch{5}.spm.spatial.preproc.tissue(4).tpm = {['C:\Program Files\spm12\tpm\TPM.nii,4']};
   matlabbatch{5}.spm.spatial.preproc.tissue(4).ngaus = 3;
   matlabbatch{5}.spm.spatial.preproc.tissue(4).native = [1 0];
   matlabbatch{5}.spm.spatial.preproc.tissue(4).warped = [0 0];
   % Soft Tissue (from SPM)
   matlabbatch{5}.spm.spatial.preproc.tissue(5).tpm = {['C:\Program Files\spm12\tpm\TPM.nii,5']};
   matlabbatch{5}.spm.spatial.preproc.tissue(5).ngaus = 4;
   matlabbatch{5}.spm.spatial.preproc.tissue(5).native = [1 0];
   matlabbatch{5}.spm.spatial.preproc.tissue(5).warped = [0 0];
   % Other Things (from SPM)
   matlabbatch{5}.spm.spatial.preproc.tissue(6).tpm = {['C:\Program Files\spm12\tpm\TPM.nii,6']};
   matlabbatch{5}.spm.spatial.preproc.tissue(6).ngaus = 2;
   matlabbatch{5}.spm.spatial.preproc.tissue(6).native = [0 0];
   matlabbatch{5}.spm.spatial.preproc.tissue(6).warped = [0 0];
   matlabbatch{5}.spm.spatial.preproc.warp.mrf = 1;
   matlabbatch{5}.spm.spatial.preproc.warp.cleanup = 1;
   matlabbatch{5}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
   matlabbatch{5}.spm.spatial.preproc.warp.affreg = 'mni';
   matlabbatch{5}.spm.spatial.preproc.warp.fwhm = 0;
   matlabbatch{5}.spm.spatial.preproc.warp.samp = 3;
   matlabbatch{5}.spm.spatial.preproc.warp.write = [0 1];

   % MODULE 6: Normalize (Write)
   matlabbatch{6}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
   % Selects deformation field in anat from Segmenting step; y_rsub.nii
   matlabbatch{6}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
   matlabbatch{6}.spm.spatial.normalise.write.subj.resample(2) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 2)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{2}, '.','files'));
   % Selects slice-time corrected images (module 2) session 1-2; arsub.nii
   matlabbatch{6}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                             78 76 85];
   matlabbatch{6}.spm.spatial.normalise.write.woptions.vox = [3 3 3];
   % voxel res of warped images, 2x2x2 is default but 3x3x3 is lower\uses less space
   matlabbatch{6}.spm.spatial.normalise.write.woptions.interp = 4;
   % Interpolates deformation field\func images via 4th degree B-Spline
   matlabbatch{6}.spm.spatial.normalise.write.woptions.prefix = 'w';
    % MODULE 7: Smoothing
   matlabbatch{7}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
   % Selects normalized images from previous step; warsub...bold.nii
   matlabbatch{7}.spm.spatial.smooth.fwhm = [8 8 8];
   % FWHM smoothin kernel in mm
   matlabbatch{7}.spm.spatial.smooth.dtype = 0;
   % Keeps output images as same datatype as original
   matlabbatch{7}.spm.spatial.smooth.im = 0;
   % Implicit masking set to 'no'
   matlabbatch{7}.spm.spatial.smooth.prefix = 's';
    % MODULE 8: File Set Split
   matlabbatch{8}.cfg_basicio.file_dir.file_ops.cfg_file_split.name = 'run1run2GramFileSplit';
   matlabbatch{8}.cfg_basicio.file_dir.file_ops.cfg_file_split.files(1) = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
   % Selects smoothed images from previous module; swarsub.nii
   matlabbatch{8}.cfg_basicio.file_dir.file_ops.cfg_file_split.index = {
                                                                        1
                                                                        2
                                                                        }';
   % Splits smoothed images into runs 1 and 2
   %MODULE 9: Model Specification
   matlabbatch{9}.spm.stats.fmri_spec.dir = {['D:\ds003604\',subject,'\ses-7\1stLevelG']};
   matlabbatch{9}.spm.stats.fmri_spec.timing.units = 'secs';
   matlabbatch{9}.spm.stats.fmri_spec.timing.RT = 1.25;
   matlabbatch{9}.spm.stats.fmri_spec.timing.fmri_t = 16;
   matlabbatch{9}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
  
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).scans(1) = cfg_dep('File Set Split: run1run2GramFileSplit (1)', substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('{}',{1}));
   % Selects Session 1 files from File Set Split Module
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).cond(1).name = 'G_Ctrl';
   data_g_ctrl_run1 = load(['D:\ds003604\' subject '\ses-7\func\timing\',subject,'_ses-7_G_Ctrl.txt']);
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).cond(1).onset = data_g_ctrl_run1(:,1);
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).cond(1).duration = data_g_ctrl_run1(:,2);
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).cond(1).tmod = 0;
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).cond(1).orth = 1;
  
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).cond(2).name = 'G_F';
   data_g_f_run1 = load(['D:\ds003604\' subject '\ses-7\func\timing\',subject,'_ses-7_G_F.txt']);
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).cond(2).onset = data_g_f_run1(:,1);
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).cond(2).duration = data_g_f_run1(:,2);
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).cond(2).tmod = 0;
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).cond(2).orth = 1;
  
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).cond(3).name = 'G_G';
   data_g_g_run1 = load(['D:\ds003604\' subject '\ses-7\func\timing\',subject,'_ses-7_G_G.txt']);
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).cond(3).onset = data_g_g_run1(:,1);
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).cond(3).duration = data_g_g_run1(:,2);
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).cond(3).tmod = 0;
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).cond(3).orth = 1;
  
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).cond(4).name = 'G_P';
   data_g_p_run1 = load(['D:\ds003604\' subject '\ses-7\func\timing\',subject,'_ses-7_G_P.txt']);
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).cond(4).onset = data_g_p_run1(:,1);
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).cond(4).duration = data_g_p_run1(:,2);
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).cond(4).tmod = 0;
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).cond(4).orth = 1;
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).multi = {''};
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).multi_reg = {''};
   matlabbatch{9}.spm.stats.fmri_spec.sess(1).hpf = 128;
  
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).scans(1) = cfg_dep('File Set Split: run1run2GramFileSplit (2)', substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('{}',{2}));
   % Selects Session 2 files from File Set Split Module
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).cond(1).name = 'G_Ctrl';
   data_g_ctrl_run2 = load(['D:\ds003604\' subject '\ses-7\func\timing\',subject,'_ses-7_G_Ctrl_2.txt']);
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).cond(1).onset = data_g_ctrl_run2(:,1);
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).cond(1).duration = data_g_ctrl_run2(:,2);
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).cond(1).tmod = 0;
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).cond(1).orth = 1;
  
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).cond(2).name = 'G_F';
   data_g_f_run2 = load(['D:\ds003604\' subject '\ses-7\func\timing\',subject,'_ses-7_G_F_2.txt']);
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).cond(2).onset = data_g_f_run2(:,1);
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).cond(2).duration = data_g_f_run2(:,2);
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).cond(2).tmod = 0;
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).cond(2).orth = 1;
  
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).cond(3).name = 'G_G';
   data_g_g_run2 = load(['D:\ds003604\' subject '\ses-7\func\timing\',subject,'_ses-7_G_G_2.txt']);
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).cond(3).onset = data_g_g_run2(:,1);
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).cond(3).duration = data_g_g_run2(:,2);
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).cond(3).tmod = 0;
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).cond(3).orth = 1;
  
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).cond(4).name = 'G_P';
   data_g_p_run2 = load(['D:\ds003604\' subject '\ses-7\func\timing\',subject,'_ses-7_G_P_2.txt']);
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).cond(4).onset = data_g_p_run2(:,1);
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).cond(4).duration = data_g_p_run2(:,2);
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).cond(4).tmod = 0;
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).cond(4).orth = 1;
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).multi = {''};
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).regress = struct('name', {}, 'val', {});
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).multi_reg = {''};
   matlabbatch{9}.spm.stats.fmri_spec.sess(2).hpf = 128;
  
   matlabbatch{9}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
   matlabbatch{9}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
   matlabbatch{9}.spm.stats.fmri_spec.volt = 1;
   matlabbatch{9}.spm.stats.fmri_spec.global = 'None';
   matlabbatch{9}.spm.stats.fmri_spec.mthresh = 0.8;
   matlabbatch{9}.spm.stats.fmri_spec.mask = {''};
   matlabbatch{9}.spm.stats.fmri_spec.cvi = 'AR(1)';
   % MODULE 10: Model Estimation
   matlabbatch{10}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
   % Selects SPM.mat file from previous module
   matlabbatch{10}.spm.stats.fmri_est.write_residuals = 0;
   % Writes residuals
   matlabbatch{10}.spm.stats.fmri_est.method.Classical = 1;
   % MODULE 11: Contrast Manager
   matlabbatch{11}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{10}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
   % Selects SPM.mat file from model estimation step
   matlabbatch{11}.spm.stats.con.consess{1}.tcon.name = 'G_F-Ctrl';
   matlabbatch{11}.spm.stats.con.consess{1}.tcon.weights = [-1 1 0 0];
   matlabbatch{11}.spm.stats.con.consess{1}.tcon.sessrep = 'replsc';
   % Repeats contrast for both sessions and weights accordingly
  
   matlabbatch{11}.spm.stats.con.consess{2}.tcon.name = 'G_G-Ctrl';
   matlabbatch{11}.spm.stats.con.consess{2}.tcon.weights = [-1 0 1 0];
   matlabbatch{11}.spm.stats.con.consess{2}.tcon.sessrep = 'replsc';
  
   matlabbatch{11}.spm.stats.con.consess{3}.tcon.name = 'G_P-Ctrl';
   matlabbatch{11}.spm.stats.con.consess{3}.tcon.weights = [-1 0 0 1];
   matlabbatch{11}.spm.stats.con.consess{3}.tcon.sessrep = 'replsc';
  
   matlabbatch{11}.spm.stats.con.consess{4}.tcon.name = 'G_F';
   matlabbatch{11}.spm.stats.con.consess{4}.tcon.weights = [0 1 0 0];
   matlabbatch{11}.spm.stats.con.consess{4}.tcon.sessrep = 'replsc';
  
   matlabbatch{11}.spm.stats.con.consess{5}.tcon.name = 'G_G';
   matlabbatch{11}.spm.stats.con.consess{5}.tcon.weights = [0 0 1 0];
   matlabbatch{11}.spm.stats.con.consess{5}.tcon.sessrep = 'replsc';
  
   matlabbatch{11}.spm.stats.con.consess{6}.tcon.name = 'G_P';
   matlabbatch{11}.spm.stats.con.consess{6}.tcon.weights = [0 0 0 1];
   matlabbatch{11}.spm.stats.con.consess{6}.tcon.sessrep = 'replsc';
   matlabbatch{11}.spm.stats.con.delete = 0;
  
  

