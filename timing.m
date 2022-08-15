% % Converts timing files from BIDS format into a two-column format that can
% % be read by SPM
% 
% % The columns are:
% % 1. Onset (in seconds); and
% % 2. Duration (in seconds


clc;
clear all;
rootdir="F:\ds003604";
%Set as directory with fMRI subject folders
d=dir(rootdir); 
cd(rootdir);
% cd('C:\Program Files\spm12');
%spm fmri

for i=209:246
    %16-247
    
    subject = num2str(d(i).name); % Zero-pads each number so that the subject ID is 2 characters long
    filepath = [subject, '/ses-7/func']
     if exist(filepath, 'dir') == 0
       continue;
     else
     cd(filepath) % Navigate to the subject's directory

    % Specify the folder where the files live.

    sem_file= ls([subject '_ses-7_task-Sem_*_run-01_events.tsv']);
    sem_onsetTimes = tdfread(sem_file, '\t') % Read onset times file
    sem_onsetTimes.trial_type = string(sem_onsetTimes.trial_type); % Convert char array to string array, to make logical comparisons easier

    S_H = [];
    S_L = [];
    S_U = [];
    S_Ctrl = [];

    for i = 1:length(sem_onsetTimes.onset)
        if strtrim(sem_onsetTimes.trial_type(i,:)) == 'S_U'
            S_U = [S_U; sem_onsetTimes.onset(i,:) sem_onsetTimes.duration(i,:)];
        elseif strtrim(sem_onsetTimes.trial_type(i,:)) == 'S_L'
            S_L = [S_L; sem_onsetTimes.onset(i,:) sem_onsetTimes.duration(i,:)];
        elseif strtrim(sem_onsetTimes.trial_type(i,:)) == 'S_H'
            S_H = [S_H; sem_onsetTimes.onset(i,:) sem_onsetTimes.duration(i,:)];
        elseif strtrim(sem_onsetTimes.trial_type(i,:)) == 'S_C'
            S_Ctrl = [S_Ctrl; sem_onsetTimes.onset(i,:) sem_onsetTimes.duration(i,:)];
        end
    end


    plaus_file= ls([subject '_ses-7_task-Plaus_*_run-01_events.tsv']);
    plaus_onsetTimes = tdfread(plaus_file, '\t') % Read onset times file
    plaus_onsetTimes.trial_type = string(plaus_onsetTimes.trial_type); % Convert char array to string array, to make logical comparisons easier

    P_Ctrl = [];
    P_I = [];
    P_S = [];
    P_W = [];

    for i = 1:length(plaus_onsetTimes.onset)
        if strtrim(plaus_onsetTimes.trial_type(i,:)) == 'SP_C'
            P_Ctrl = [P_Ctrl; plaus_onsetTimes.onset(i,:) plaus_onsetTimes.duration(i,:)];
        elseif strtrim(plaus_onsetTimes.trial_type(i,:)) == 'SP_I'
            P_I = [P_I; plaus_onsetTimes.onset(i,:) plaus_onsetTimes.duration(i,:)];
        elseif strtrim(plaus_onsetTimes.trial_type(i,:)) == 'SP_S'
            P_S = [P_S; plaus_onsetTimes.onset(i,:) plaus_onsetTimes.duration(i,:)];
        elseif strtrim(plaus_onsetTimes.trial_type(i,:)) == 'SP_W'
            P_W = [P_W; plaus_onsetTimes.onset(i,:) plaus_onsetTimes.duration(i,:)];
        end
    end

    gram_file= ls([subject '_ses-7_task-Gram_*_run-01_events.tsv']);
    gram_onsetTimes = tdfread(gram_file, '\t') % Read onset times file
    gram_onsetTimes.trial_type = string(gram_onsetTimes.trial_type); % Convert char array to string array, to make logical comparisons easier

    G_Ctrl = [];
    G_F = [];
    G_G = [];
    G_P = [];

    for i = 1:length(gram_onsetTimes.onset)
        if strtrim(gram_onsetTimes.trial_type(i,:)) == 'G_C'
            G_Ctrl = [G_Ctrl; gram_onsetTimes.onset(i,:) gram_onsetTimes.duration(i,:)];
        elseif strtrim(gram_onsetTimes.trial_type(i,:)) == 'G_F'
            G_F = [G_F; gram_onsetTimes.onset(i,:) gram_onsetTimes.duration(i,:)];
        elseif strtrim(gram_onsetTimes.trial_type(i,:)) == 'G_G'
            G_G = [G_G; gram_onsetTimes.onset(i,:) gram_onsetTimes.duration(i,:)];
        elseif strtrim(gram_onsetTimes.trial_type(i,:)) == 'G_P'
            G_P = [G_P; gram_onsetTimes.onset(i,:) gram_onsetTimes.duration(i,:)];
        end
    end

    % Save timing files into text files

    save([subject,'_ses-7_S_Ctrl.txt'], 'S_Ctrl', '-ASCII');
    save([subject,'_ses-7_S_H.txt'], 'S_H', '-ASCII');
    save([subject,'_ses-7_S_U.txt'], 'S_U', '-ASCII');
    save([subject,'_ses-7_S_L.txt'], 'S_L', '-ASCII');

    save([subject,'_ses-7_P_Ctrl.txt'], 'P_Ctrl', '-ASCII');
    save([subject,'_ses-7_P_I.txt'], 'P_I', '-ASCII');
    save([subject,'_ses-7_P_S.txt'], 'P_S', '-ASCII');
    save([subject,'_ses-7_P_W.txt'], 'P_W', '-ASCII');

    save([subject,'_ses-7_G_Ctrl.txt'], 'G_Ctrl', '-ASCII');
    save([subject,'_ses-7_G_F.txt'], 'G_F', '-ASCII');
    save([subject,'_ses-7_G_G.txt'], 'G_G', '-ASCII');
    save([subject,'_ses-7_G_P.txt'], 'G_P', '-ASCII');
   
    % Go back to Subject directory
     end

    clear G_* gram* P_* plaus* S_* sem*
    cd ../../..

end

%   [rs,cs]=size(sem_files)
%         min_num=1000;
%         smin_file=[];
%         num=0;
%         for file=1:rs
%             D=str2double(extractBetween(sem_files(file,:),"D","S"))
%             S=str2double(extractBetween(sem_files(file,:),['D',num2str(D),'S'],"_run"))
%             num=D+S;
%             if num<min_num==0
%                 continue;
%             else
%                 smin_file=sem_files(file,:);
%                 min_num=num;
%             end
%             clear D S num
%         end



