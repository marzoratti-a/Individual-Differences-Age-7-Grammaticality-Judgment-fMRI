clc;
clear all;
rootdir="F:\ds003604\";
%Set as directory with fMRI subject folders
d=dir(rootdir); 
cd(rootdir);
% cd('C:\Program Files\spm12');
%spm fmri

% clc;
% clear all;
% rootdir="C:\Users\anama\OneDrive\Documents\LabWork\proj-62bddf5ef3194eded6f9293d\bids\derivatives\upload\";
% %Set as directory with fMRI subject folders
% d=dir(rootdir); 
% cd(rootdir);
% 
% for i= 18
%      filepath = [d(i).name, '\ses-7\func']
%      cd(filepath);
%      %delete('*_Head_tmp.txt')
%      %delete('*_all_tmp.txt')
%      gunzip('*.gz')
%      cd(rootdir);
%  end 

 for i= 17:246
     %func folders up to sub-5078/i=66 unzipped for ses-5
    filepath = [char(d(i).name), '\ses-7\func'];
    if exist(filepath, 'dir') == 0
       continue;
    else
       cd(filepath); 
     % delete([char(d(i).name),'_ses-9_task-Phon*'])
%        gunzip('*.gz')
      time= [filepath,'\timing']
       if exist(time, 'dir')==0
           mkdir timing
       else
           continue;
       end
       movefile('*.txt','timing')
    end
    cd(rootdir);

end 

