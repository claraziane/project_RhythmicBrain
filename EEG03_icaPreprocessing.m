%% Preprocessing: phase 3
% - Run ICA (once for all conditions)
% - Save with ICA weights

close all;
clear all;
clc;

% Declare paths
pathData = '/Volumes/Seagate/project_rhythmicBrain/DATA/';
addpath('/Users/claraziane/Documents/Académique/Informatique/MATLAB/eeglab2021.1');

Participants   = {'sub-001'; 'sub-002'; 'sub-003'; 'sub-004'; 'sub-005'; 'sub-006'; 'sub-007'; 'sub-008'; 'sub-009'; 'sub-010'; 'sub-011'; 'sub-012'; 'sub-013'; 'sub-014'; 'sub-015'; 'sub-016'; 'sub-017'; 'sub-018'};
Blocks         = {'run-01'; 'run-02'; 'run-03'; 'run-04'; 'run-05'; 'run-06'; 'run-07'; 'run-08'; 'run-09'; 'run-10'; 'run-11'; 'run-12'; 'run-13'; 'run-14'; 'run-15'; 'run-16'; 'run-17'; 'run-18'; 'run-19'};
nBlocks        = 0;

extensionRoot  = '_chanClean.set';
extensionFinal = '_ica.set';

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
for iParticipant = 1:length(Participants)
    disp(Participants{iParticipant});
    directory = fullfile(pathData, Participants{iParticipant}, 'eeg');

    tic
    for iBlock = 1:length(Blocks)
        fileRead = [Participants{iParticipant} '_task-AudioCueWalkingStudy_' Blocks{iBlock} extensionRoot];
        
        % To deal with different snumber of blocks per participant
        if exist(fullfile(directory, fileRead), 'file') ~= 0 
            EEG = pop_loadset('filename', fileRead,'filepath', directory); % Loads an EEG dataset
            [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, 0); % Edits/saves EEG dataset structure information
            
            nBlocks = nBlocks + 1;
        end

    end

    % Run ICA on all conditions at once
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, nBlocks,'retrieve',[1:nBlocks] ,'study',0);
    EEG = eeg_checkset(EEG);
    EEG = pop_runica(EEG, 'icatype','runica','concatcond','on', 'chanind',[1:108], 'options',{'extended',1});
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset(EEG);

    % Save each condition pruned with ICA as separate file
    for iBlock = 1:nBlocks
        fileWrite    = [Participants{iParticipant} '_task-AudioCueWalkingStudy_' Blocks{iBlock} extensionFinal];

        % Retrieve condition to save
        if iBlock == 1
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, [iBlock:nBlocks] ,'retrieve', iBlock ,'study',0);
        else
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, iBlock-1 ,'retrieve', iBlock ,'study',0);

        end
        EEG = eeg_checkset(EEG);
        toc

        % Save
        EEG = pop_saveset(EEG, 'filename', fileWrite, 'filepath', [directory '/']); %pop_saveset saves one or more EEG dataset structures
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

    end

    ALLEEG = [];

end