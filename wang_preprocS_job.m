  subject = num2str(d(i).name,'%02d'); % Zero-pads each number so that the subject ID is 2 characters long

    %%%%%%%%%%
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.name = [subject,'_Semruns'];
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {
                                                                        {['D:\ds003604\',subject,'\ses-7\func\',subject,'_ses-7_task-Sem_acq_run-01_bold.nii']}
                                                                        {['D:\ds003604\',subject,'\ses-7\func\',subject,'_ses-7_task-Sem_acq_run-02_bold.nii']}
                                                                        }';
    % MODULE 2: Realign (Estimate & Reslice)
    matlabbatch{2}.spm.spatial.realign.estwrite.data{1}(1) = cfg_dep(['Named File Selector: ',subject,'_Semruns(1) - Files'], substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
    matlabbatch{2}.spm.spatial.realign.estwrite.data{2}(1) = cfg_dep(['Named File Selector: ',subject,'_Semruns(2) - Files'], substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{2}));
    % Specifies run 1 and 2 files selected in named file selector module ^^^
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.fwhm = 5; 
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

  % MODULE 3: Segmenting
    matlabbatch{3}.spm.spatial.preproc.channel.vols = {['D:\ds003604\',subject,'\ses-7\anat_S\',subject,'_ses-7_acq_T1w.nii,1']};
    matlabbatch{3}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{3}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{3}.spm.spatial.preproc.channel.write = [0 1];
    matlabbatch{3}.spm.spatial.preproc.tissue(1).tpm = {'C:\Users\anama\OneDrive\Documents\My Experiments\demos\mw_com_prior_Age_0084.nii,1'};
    matlabbatch{3}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{3}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{3}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{3}.spm.spatial.preproc.tissue(2).tpm = {'C:\Users\anama\OneDrive\Documents\My Experiments\demos\mw_com_prior_Age_0084.nii,2'};
    matlabbatch{3}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{3}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{3}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{3}.spm.spatial.preproc.tissue(3).tpm = {'C:\Users\anama\OneDrive\Documents\My Experiments\demos\mw_com_prior_Age_0084.nii,3'};
    matlabbatch{3}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{3}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{3}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{3}.spm.spatial.preproc.tissue(4).tpm = {'C:\Users\anama\OneDrive\Documents\My Experiments\demos\mw_com_prior_Age_0084.nii,4'};
    matlabbatch{3}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{3}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{3}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{3}.spm.spatial.preproc.tissue(5).tpm = {'C:\Users\anama\OneDrive\Documents\My Experiments\demos\mw_com_prior_Age_0084.nii,5'};
    matlabbatch{3}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{3}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{3}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{3}.spm.spatial.preproc.tissue(6).tpm = {'C:\Users\anama\OneDrive\Documents\My Experiments\demos\mw_com_prior_Age_0084.nii,6'};
    matlabbatch{3}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{3}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{3}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{3}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{3}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{3}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{3}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{3}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{3}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{3}.spm.spatial.preproc.warp.write = [0 1];
    matlabbatch{3}.spm.spatial.preproc.warp.vox = NaN;
    matlabbatch{3}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
                                                  NaN NaN NaN];
    
    % MODULE 4: Skull Stripping
    matlabbatch{4}.spm.util.imcalc.input(1) = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
    matlabbatch{4}.spm.util.imcalc.input(2) = cfg_dep('Segment: c2 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
    matlabbatch{4}.spm.util.imcalc.input(3) = cfg_dep('Segment: c3 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));
    matlabbatch{4}.spm.util.imcalc.output = 'mask';,
    matlabbatch{4}.spm.util.imcalc.outdir = {['D:\ds003604\',subject,'\ses-7\anat_S']};
    matlabbatch{4}.spm.util.imcalc.expression = 'i1+i2+i3';
    matlabbatch{4}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{4}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{4}.spm.util.imcalc.options.mask = 0;
    matlabbatch{4}.spm.util.imcalc.options.interp = 0;
    matlabbatch{4}.spm.util.imcalc.options.dtype = 16;
    
    % MODULE 5: Mask Making
    matlabbatch{5}.spm.util.imcalc.input(1) = cfg_dep('Segment: Bias Corrected (1)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
    matlabbatch{5}.spm.util.imcalc.input(2) = cfg_dep('Image Calculator: ImCalc Computed Image: mask', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{5}.spm.util.imcalc.output = 'noskull';
    matlabbatch{5}.spm.util.imcalc.outdir = {['D:\ds003604\',subject,'\ses-7\anat_S']};
    matlabbatch{5}.spm.util.imcalc.expression = 'i1.*(i2>0)';
    matlabbatch{5}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{5}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{5}.spm.util.imcalc.options.mask = 0;
    matlabbatch{5}.spm.util.imcalc.options.interp = 0;
    matlabbatch{5}.spm.util.imcalc.options.dtype = 16;
    

  % MODULE 6: Coregister (Estimate & Reslice)
    matlabbatch{6}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
    matlabbatch{6}.spm.spatial.coreg.estimate.source(1) = cfg_dep('Image Calculator: ImCalc Computed Image: noskull', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{6}.spm.spatial.coreg.estimate.other = {''};
    matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

   % MODULE 7: Normalize
    matlabbatch{7}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
    matlabbatch{7}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rfiles'));
    matlabbatch{7}.spm.spatial.normalise.write.subj.resample(2) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 2)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{2}, '.','rfiles'));
    matlabbatch{7}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                              78 76 85];
    matlabbatch{7}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
    matlabbatch{7}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{7}.spm.spatial.normalise.write.woptions.prefix = 'w';

  % MODULE 8: Smooth
    matlabbatch{8}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    matlabbatch{8}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{8}.spm.spatial.smooth.dtype = 0;
    matlabbatch{8}.spm.spatial.smooth.im = 0;
    matlabbatch{8}.spm.spatial.smooth.prefix = 's';

    % MODULE 9: File Set Split
    matlabbatch{9}.cfg_basicio.file_dir.file_ops.cfg_file_split.name = 'run1run2SemFileSplit';
    matlabbatch{9}.cfg_basicio.file_dir.file_ops.cfg_file_split.files(1) = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    % Selects smoothed images from previous module; swarsub.nii
    matlabbatch{9}.cfg_basicio.file_dir.file_ops.cfg_file_split.index = {
                                                                         1
                                                                         2
                                                                         }';
    % Splits smoothed images into runs 1 and 2
    %MODULE 10: Model Specification
    matlabbatch{10}.spm.stats.fmri_spec.dir = {['D:\ds003604\',subject,'\ses-7\1stLevelS']};
    matlabbatch{10}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{10}.spm.stats.fmri_spec.timing.RT = 1.25;
    matlabbatch{10}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{10}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).scans(1) = cfg_dep('File Set Split: run1run2SemFileSplit (1)', substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('{}',{1}));
    % Selects Session 1 files from File Set Split Module
    % Conditions
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).cond(1).name = 'S_Ctrl';
    data_s_ctrl_run1 = load(['D:\ds003604\',subject,'\ses-7\func\timing\',subject,'_ses-7_S_Ctrl.txt']);
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).cond(1).onset = data_s_ctrl_run1(:,1);
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).cond(1).duration = data_s_ctrl_run1(:,2);
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).cond(1).tmod = 0;
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).cond(1).orth = 1;
    
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).cond(2).name = 'S_H';
    data_s_h_run1 = load(['D:\ds003604\' subject '\ses-7\func\timing\',subject,'_ses-7_S_H.txt']);
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).cond(2).onset = data_s_h_run1(:,1);
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).cond(2).duration = data_s_h_run1(:,2);
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).cond(2).tmod = 0;
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).cond(2).orth = 1;
    
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).cond(3).name = 'S_L';
    data_s_l_run1 = load(['D:\ds003604\' subject '\ses-7\func\timing\',subject,'_ses-7_S_L.txt']);
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).cond(3).onset = data_s_l_run1(:,1);
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).cond(3).duration = data_s_l_run1(:,2);
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).cond(3).tmod = 0;
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).cond(3).orth = 1;
    
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).cond(4).name = 'S_U';
    data_s_u_run1 = load(['D:\ds003604\' subject '\ses-7\func\timing\',subject,'_ses-7_S_U.txt']);
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).cond(4).onset = data_s_u_run1(:,1);
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).cond(4).duration = data_s_u_run1(:,2);
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).cond(4).tmod = 0;
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).cond(4).orth = 1;
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).multi = {''};
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).multi_reg = {''};
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).hpf = 128;
    
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).scans(1) = cfg_dep('File Set Split: run1run2SemFileSplit (2)', substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('{}',{2}));
    % Selects Session 2 files from File Set Split Module
    % Conditions
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).cond(1).name = 'S_Ctrl';
    data_s_ctrl_run2 = load(['D:\ds003604\' subject '\ses-7\func\timing\',subject,'_ses-7_S_Ctrl_2.txt']);
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).cond(1).onset = data_s_ctrl_run2(:,1);
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).cond(1).duration = data_s_ctrl_run2(:,2);
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).cond(1).tmod = 0;
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).cond(1).orth = 1;
    
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).cond(2).name = 'S_H';
    data_s_h_run2 = load(['D:\ds003604\' subject '\ses-7\func\timing\',subject,'_ses-7_S_H_2.txt']);
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).cond(2).onset = data_s_h_run2(:,1);
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).cond(2).duration = data_s_h_run2(:,2);
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).cond(2).tmod = 0;
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).cond(2).orth = 1;
    
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).cond(3).name = 'S_L';
    data_s_l_run2 = load(['D:\ds003604\' subject '\ses-7\func\timing\',subject,'_ses-7_S_L_2.txt']);
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).cond(3).onset = data_s_l_run2(:,1);
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).cond(3).duration = data_s_l_run2(:,2);
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).cond(3).tmod = 0;
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).cond(3).orth = 1;
    
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).cond(4).name = 'S_U';
    data_s_u_run2 = load(['D:\ds003604\' subject '\ses-7\func\timing\',subject,'_ses-7_S_U_2.txt']);
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).cond(4).onset = data_s_u_run2(:,1);
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).cond(4).duration = data_s_u_run2(:,2);
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).cond(4).tmod = 0;
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).cond(4).orth = 1;
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).multi = {''};
    matlabbatch{10}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).multi_reg = {''};
    matlabbatch{10}.spm.stats.fmri_spec.sess(2).hpf = 128;
    matlabbatch{10}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{10}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{10}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{10}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{10}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{10}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{10}.spm.stats.fmri_spec.cvi = 'AR(1)';

    % MODULE 11: Model Estimation
    matlabbatch{11}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{10}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    % Selects SPM.mat file from previous module
    matlabbatch{11}.spm.stats.fmri_est.write_residuals = 0;
    % Writes residuals
    matlabbatch{11}.spm.stats.fmri_est.method.Classical = 1;

    % MODULE 12: Contrast Manager
    matlabbatch{12}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{11}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    % Selects SPM.mat file from model estimation step
    % Selects SPM.mat file from model estimation step
    matlabbatch{12}.spm.stats.con.consess{1}.tcon.name = 'S_H-Ctrl';
    matlabbatch{12}.spm.stats.con.consess{1}.tcon.weights = [-1 1 0 0];
    matlabbatch{12}.spm.stats.con.consess{1}.tcon.sessrep = 'replsc';
    % Repeats contrast for both sessions and weights accordingly
    matlabbatch{12}.spm.stats.con.consess{2}.tcon.name = 'S_L-Ctrl';
    matlabbatch{12}.spm.stats.con.consess{2}.tcon.weights = [-1 0 1 0];
    matlabbatch{12}.spm.stats.con.consess{2}.tcon.sessrep = 'replsc';
    matlabbatch{12}.spm.stats.con.consess{3}.tcon.name = 'S_U-Ctrl';
    matlabbatch{12}.spm.stats.con.consess{3}.tcon.weights = [-1 0 0 1];
    matlabbatch{12}.spm.stats.con.consess{3}.tcon.sessrep = 'replsc';
    matlabbatch{12}.spm.stats.con.consess{4}.tcon.name = 'S_H';
    matlabbatch{12}.spm.stats.con.consess{4}.tcon.weights = [0 1 0 0];
    matlabbatch{12}.spm.stats.con.consess{4}.tcon.sessrep = 'replsc';
    matlabbatch{12}.spm.stats.con.consess{5}.tcon.name = 'S_L';
    matlabbatch{12}.spm.stats.con.consess{5}.tcon.weights = [0 0 1 0];
    matlabbatch{12}.spm.stats.con.consess{5}.tcon.sessrep = 'replsc';
    matlabbatch{12}.spm.stats.con.consess{6}.tcon.name = 'S_U';
    matlabbatch{12}.spm.stats.con.consess{6}.tcon.weights = [0 0 0 1];
    matlabbatch{12}.spm.stats.con.consess{6}.tcon.sessrep = 'replsc';
    matlabbatch{12}.spm.stats.con.delete = 0;