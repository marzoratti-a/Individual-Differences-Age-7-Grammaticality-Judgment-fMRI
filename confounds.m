
% Making Regressor Files
clc;
clear all;
rootdir="C:\Users\anama\OneDrive\Documents\LabWork\proj-62bddf5ef3194eded6f9293d\bids\derivatives\upload\";
% 
%Set as directory with fMRI subject folders
d=dir(rootdir); 
cd(rootdir);
% cd('C:\Program Files\spm12');
%spm fmri

 for i=[101, 120]
    % 
%     subject = extractBefore(char(d(i).name),'.ses-7'); 
%     cd(['C:\Users\anama\OneDrive\Documents\LabWork\proj-62bddf5ef3194eded6f9293d\',subject,'.ses-7']);
%      regs= cellstr(ls(['C:\Users\anama\OneDrive\Documents\LabWork\proj-62bddf5ef3194eded6f9293d\',subject,'.ses-7\dt-neuro-regressors*']));
%      newdir=['C:\Users\anama\OneDrive\Documents\LabWork\proj-62bddf5ef3194eded6f9293d\bids\derivatives\upload\',subject,'\ses-7\func\'];
%      pat= {'sem','plaus' 'gram'};
%      for n= 1:length(regs)
%          cd(char(regs(n)));
%          fid  = fopen('_info.json','r');
%          data = textscan(fid,'%s');
%          fclose(fid);
%          Str = string(data{:});
%          if ~exist([newdir,'confounds'], 'dir')
%              mkdir([newdir,'confounds']);
%          end
%          for z= 1:3
%              if sum(double(contains(Str, pat(z))))==0
%                  continue
%              else
%                 newreg= [char(pat(z)),'_regressors.tsv'];
%                 fullfile= [newdir,'confounds\',newreg]
%                 copyfile ('regressors.tsv', fullfile);
%              end
%              clear newreg fullfile
%          end 
%          cd ..
%      end 
%  end

     subject= char(d(i).name)
     cd([subject,'\ses-7\func\confounds']);
     task={'sem','plaus','gram'}
     for n=1:3
         csf= importdata(['csf_',char(task(n)),'_noHead_tmp.txt']);
         csf(:,1);
         global_signal= importdata(['global_signal_',char(task(n)),'_noHead_tmp.txt']);
         global_signal(:,1);
         rot_x= importdata(['rot_x_',char(task(n)),'_noHead_tmp.txt']);
         rot_x(:,1);
         rot_y= importdata(['rot_y_',char(task(n)),'_noHead_tmp.txt']);
         rot_y(:,1);
         rot_z= importdata(['rot_z_',char(task(n)),'_noHead_tmp.txt']);
         rot_z(:,1);
         trans_x= importdata(['trans_x_',char(task(n)),'_noHead_tmp.txt']);
         trans_x(:,1);
         trans_y= importdata(['trans_y_',char(task(n)),'_noHead_tmp.txt']);
         trans_y(:,1);
         trans_z= importdata(['trans_z_',char(task(n)),'_noHead_tmp.txt']);
         trans_z(:,1);
         white_matter= importdata(['white_matter_',char(task(n)),'_noHead_tmp.txt']);
         white_matter(:,1);
     
         names=cell(9,1);
            names(1)= {'csf'};
            names(2)= {'global_signal'};
            names(3)= {'rot_x'};
            names(4)= {'rot_y'};
            names(5)= {'rot_z'};
            names(6)= {'trans_x'};
            names(7)= {'trans_y'};
            names(8)= {'trans_z'};
            names(9)= {'white_matter'};
            
         regressors=cell(9,1);
            regressors(1)= {csf};
            regressors(2)= {global_signal};
            regressors(3)= {rot_x};
            regressors(4)= {rot_y};
            regressors(5)= {rot_z};
            regressors(6)= {trans_x};
            regressors(7)= {trans_y};
            regressors(8)= {trans_z};
            regressors(9)= {white_matter};

            reg_matrix= [names, regressors];
            save([char(task(n)),'_reg.mat'], 'reg_matrix')

            clear reg_matrix csf global_signal rot_x rot_y rot_z trans_x trans_y trans_z white_matter;
     end
     cd ../../../..
end