clc;
clear all;
rootdir="F:\ds003604";
%Set as directory with fMRI subject folders
d=dir(rootdir); 
cd(rootdir);
%spm fmri

for i=16
    
    subject = num2str(d(i).name); % Zero-pads each number so that the subject ID is 2 characters long
    cd([subject, '/ses-7/func']) % Navigate to the subject's directory

    % Making Multiple Condition Files

    S_Ctrl = importdata([subject, '_ses-7_S_Ctrl.txt']);
    S_Ctrl(:,1);
    S_L = importdata([subject,'_ses-7_S_L.txt']);
    S_L(:,1);
    S_H = importdata([subject,'_ses-7_S_H.txt']);
    S_H(:,1);
    S_U = importdata([subject,'_ses-7_S_U.txt']);
    S_U(:,1);

    name= cell(4,1); 
        name(1)= {'S_Ctrl'};
        name(2)= {'S_L'};
        name(3)= {'S_H'};
        name(4)= {'S_U'};
    onset= cell(4,1);
        onset(1)= {S_Ctrl(:,1)};
        onset(2)= {S_L(:,1)};
        onset(3)= {S_H(:,1)};
        onset(4)= {S_U(:,1)};
    duration= cell(4,1);
        duration(1)= {S_Ctrl(:,2)};
        duration(2)= {S_L(:,2)};
        duration(3)= {S_H(:,2)};
        duration(4)= {S_U(:,2)};

    sem_matrix= [name onset duration];
    save Sem_r1.mat sem_matrix

    clear name onset duration;

    P_Ctrl = importdata([subject,'_ses-7_P_Ctrl.txt']);
    P_Ctrl(:,1);
    P_I = importdata([subject,'_ses-7_P_I.txt']);
    P_I(:,1);
    P_S = importdata([subject,'_ses-7_P_S.txt']);
    P_S(:,1);
    P_W = importdata([subject,'_ses-7_P_W.txt']);
    P_W(:,1);

    name= cell(4,1); 
        name(1)= {'P_Ctrl'};
        name(2)= {'P_I'};
        name(3)= {'P_S'};
        name(4)= {'P_W'};
    onset= cell(4,1);
        onset(1)= {P_Ctrl(:,1)};
        onset(2)= {P_I(:,1)};
        onset(3)= {P_S(:,1)};
        onset(4)= {P_W(:,1)};
    duration= cell(4,1);
        duration(1)= {P_Ctrl(:,2)};
        duration(2)= {P_I(:,2)};
        duration(3)= {P_S(:,2)};
        duration(4)= {P_W(:,2)};

    plaus_matrix= [name onset duration];
    save Plaus_r1.mat plaus_matrix

    clear name onset duration;

    G_Ctrl = importdata([subject,'_ses-7_G_Ctrl.txt']);
    G_Ctrl(:,1);
    G_F = importdata([subject,'_ses-7_G_F.txt']);
    G_F(:,1);
    G_G = importdata([subject,'_ses-7_G_G.txt']);
    G_G(:,1);
    G_P = importdata([subject,'_ses-7_G_P.txt']);
    G_P(:,1);

    name= cell(4,1); 
        name(1)= {'G_Ctrl'};
        name(2)= {'G_F'};
        name(3)= {'G_G'};
        name(4)= {'G_P'};
    onset= cell(4,1);
        onset(1)= {G_Ctrl(:,1)};
        onset(2)= {G_F(:,1)};
        onset(3)= {G_G(:,1)};
        onset(4)= {G_P(:,1)};
    duration= cell(4,1);
        duration(1)= {G_Ctrl(:,2)};
        duration(2)= {G_F(:,2)};
        duration(3)= {G_G(:,2)};
        duration(4)= {G_P(:,2)};

    gram_matrix= [name onset duration];
    save Gram_r1.mat gram_matrix

    clear name onset duration;

    cd ../../..

end