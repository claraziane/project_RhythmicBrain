%% Preprocessing : phase 1
% - Remove DC offset
% - Rereference according to the common average
% - Save
%
% Note: Data was already high-pass filtered at 1 Hz, low-pass filtered at 200 Hz
% and notch-filtered at 500 Hz before being put in open access

close all;
clear all;
clc;

% Declare paths
pathData = '/Volumes/Seagate/project_rhythmicBrain/DATA/';
addpath('/Users/claraziane/Documents/Académique/Informatique/MATLAB/eeglab2021.1')

Participants = {'sub-001'; 'sub-002'; 'sub-003'; 'sub-004'; 'sub-005'; 'sub-006'; 'sub-007'; 'sub-008'; 'sub-009'; 'sub-010'; 'sub-011'; 'sub-012'; 'sub-013'; 'sub-014'; 'sub-015'; 'sub-016'; 'sub-017'; 'sub-018'};
Blocks       = {'run-01'; 'run-02'; 'run-03'; 'run-04'; 'run-05'; 'run-06'; 'run-07'; 'run-08'; 'run-09'; 'run-10'; 'run-11'; 'run-12'; 'run-13'; 'run-14'; 'run-15'; 'run-16'; 'run-17'; 'run-18'; 'run-19'};

extensionRoot  = '_eeg.set';
extensionFinal = '.set';

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
for iParticipant = 1:length(Participants)
    directory = fullfile(pathData, Participants{iParticipant}, '/eeg/');

    for iBlock = 1:length(Blocks)

        fileRead  = [Participants{iParticipant} '_task-AudioCueWalkingStudy_' Blocks{iBlock} extensionRoot];
        fileWrite = [Participants{iParticipant} '_task-AudioCueWalkingStudy_' Blocks{iBlock} extensionFinal];

        % To deal with different number of blocks per participant
        if exist(fullfile(directory, fileRead), 'file') ~= 0 

            EEG = pop_loadset('filename', fileRead,'filepath', directory); % Loads an EEG dataset
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); % Edits/saves EEG dataset structure information
            nChan = length(EEG.chanlocs);

            % Remove baseline of the signal
            EEG = pop_rmbase(EEG, [],[]);

            % Rereferencing according to common average
            EEG = pop_reref(EEG, [], 'exclude', [109:nChan]);
            EEG = eeg_checkset(EEG);
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');

            % Save new .set file in preprocessed folder
            EEG = pop_saveset(EEG, 'filename', fileWrite, 'filepath', [directory '/']); %pop_saveset saves one or more EEG dataset structures
            ALLEEG = [];

        end

    end

end