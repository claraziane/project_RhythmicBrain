%% Preprocessing : phase 4
% - Remove ICs of occular movements
% - Save
%
% Tip : place a stopper at line 41 so that you can inspect ICs before being
% asked which ones you would like to remove

close all;
clear all;
clc;

% Declare paths
pathData = '/Volumes/Seagate/project_rhythmicBrain/DATA/';
addpath('/Users/claraziane/Documents/Académique/Informatique/MATLAB/eeglab2021.1');
load('/Volumes/Seagate/project_rhythmicBrain/Results//subAll/icaReject.mat')
load('/Volumes/Seagate/project_rhythmicBrain/DATA/chanLocs.mat') % File containing EEG electrode coordinates

Participants   = {'sub-001'; 'sub-002'; 'sub-003'; 'sub-004'; 'sub-005'; 'sub-006'; 'sub-007'; 'sub-008'; 'sub-009'; 'sub-010'; 'sub-011'; 'sub-012'; 'sub-013'; 'sub-014'; 'sub-015'; 'sub-016'; 'sub-017'; 'sub-018'};
Blocks         = {'run-01'; 'run-02'; 'run-03'; 'run-04'; 'run-05'; 'run-06'; 'run-07'; 'run-08'; 'run-09'; 'run-10'; 'run-11'; 'run-12'; 'run-13'; 'run-14'; 'run-15'; 'run-16'; 'run-17'; 'run-18'; 'run-19'};

extensionRoot  = '_ica.set';
extensionFinal = '_icaClean.set';

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
for iParticipant = 1:length(Participants)
    directory = fullfile(pathData, Participants{iParticipant}, 'eeg');
    participantStr = strcat('SUB', Participants{iParticipant}(end-2:end));

    for iBlock = 1:length(Blocks)
        conditionStr = strcat('Run', Blocks{iBlock}(end-1:end));
        fileRead     = [Participants{iParticipant} '_task-AudioCueWalkingStudy_' Blocks{iBlock} extensionRoot];
        fileWrite    = [Participants{iParticipant} '_task-AudioCueWalkingStudy_' Blocks{iBlock} extensionFinal];

        % To deal with different snumber of blocks per participant
        if exist(fullfile(directory, fileRead), 'file') ~= 0      

            EEG = pop_loadset('filename', fileRead,'filepath', directory); % Loads an EEG dataset
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','on'); % Edits/saves EEG dataset structure information
            
            iIC = 1;
            while(1)
                icaRemove = input('Should some ICAs be removed from the data ?, Y/N [Y]:', 's');
                if icaRemove == 'Y'
                    icaReject.(participantStr).(conditionStr)(iIC) = input('Which ICA should be removed from the data ?');
                    iIC = iIC + 1;
                elseif icaRemove == 'N'
                    break
                end
            end
            
            % Remove identified bad components
            if iIC > 1
                EEG = pop_subcomp(EEG, icaReject.(participantStr).(conditionStr)(:), 0); %removes ICA from EEG and subtracts their activities from the data
                eeg_checkset(EEG);
                [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');

                % Rereference
                EEG = pop_reref(EEG, []);
                [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');
            end

            % Save new .set file
            EEG = pop_saveset(EEG, 'filename', fileWrite, 'filepath', directory);
            save('/Volumes/Seagate/project_rhythmicBrain/Results/subAll/icaReject', 'icaReject');

            ALLEEG = [];

        end

    end

end