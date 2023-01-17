%% Extract relevant data from events
% Use events to extract data corresponding to five conditions : 
% (1) cued gait at preferred cadence; (2) cued gait at slow cadence; (3) cued gait at fast cadence;
% (4) uncued gait at preferred cadence; (5) standing upright while listenning to cues matching preferred cadence

clear;
close all;
clc;

% Declare paths
pathData    = '/Volumes/Seagate/project_rhythmicBrain/DATA/'; %Folder where all data is
pathResults = '/Volumes/Seagate/project_rhythmicBrain/Results/';  %Folder where to save results
addpath('/Users/claraziane/Documents/Académique/Informatique/MATLAB/eeglab2021.1'); %EEGLAB toolbox

Participants = {'sub-001'; 'sub-002'; 'sub-003'; 'sub-004'; 'sub-005'; 'sub-006'; 'sub-007'; 'sub-008'; 'sub-009'; 'sub-010'; 'sub-011'; 'sub-012'; 'sub-013'; 'sub-014'; 'sub-015'; 'sub-016'; 'sub-017'; 'sub-018'};
Blocks       = {'run-01'; 'run-02'; 'run-03'; 'run-04'; 'run-05'; 'run-06'; 'run-07'; 'run-08'; 'run-09'; 'run-10'; 'run-11'; 'run-12'; 'run-13'; 'run-14'; 'run-15'; 'run-16'; 'run-17'; 'run-18'; 'run-19'};
Conditions   = {'Pref'; 'Slow'; 'Fast'; 'uncuedPref'; 'cuedPrefRest'};
Event        = [100 200 444 999 333 114 113 119 171 271 159 259 153 253];
% "100": "Right heel strikes during uncued walking"
% "200": "Left heel strikes during uncued walking"
% "444": "First auditory cue during preferred cadence walking"
% "999": "First auditory cue marking step delay tempo shift"
% "333": "First auditory cue marking step advance tempo shift"
% "114": "First right cued step during preferred cadence walking"
% "113": "First right heel strike relative to step advance tempo shift"
% "119": "First right heel strike relative to step delay tempo shift"
% "171": "Cue during standing previously associated with right heel strike"
% "271": "Cue during standing previously associated with left heel strike"
% "159": "Auditory cues for right heel strikes during step delay pacing tempo"
% "259": "Auditory cues for left heel strikes during step delay pacing tempo"
% "153": "Auditory cues for right heel strikes during step advance pacing tempo"
% "253": "Auditory cues for left heel strikes during step advance pacing tempo"

% Electrode used for 'best-electrode' analyses (see EEG06 script)
electrode = 'cz';

% Load matrix with already-extracted data (comment if running script for the 1st time)
load([pathResults 'subAll/DATA.mat']);

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for iParticipant = 1%:length(Participants)
    participantStr = strcat('SUB', Participants{iParticipant}(end-2:end));
    k = 1; %Indexing standing (rest) trials

    for iBlock = 1:length(Blocks)

        %% Load EEG
        fileRead = ([Participants{iParticipant} '_task-AudioCueWalkingStudy_' Blocks{iBlock} '_icaClean.set']);
        directory = fullfile(pathData, Participants{iParticipant}, '/eeg/');

        % To deal with different snumber of blocks per participant
        if exist(fullfile(directory, fileRead), 'file') ~= 0
            EEG = pop_loadset('filename', fileRead, 'filepath', directory); % Loads an EEG datas
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','on'); % Edits/saves EEG dataset structure information

            % Increase data precision (from single to double) for numerical stability
            data = double(EEG.data);

            % Extract relevant info
            rateEEG.(participantStr)  = EEG.srate;
            chanLocs.(participantStr) = EEG.chanlocs;
            nChan = length(chanLocs.(participantStr));
            
            % Check that the data  matrix is full rank
            if rank(data) ~= nChan
                warning('Data matrix is not full rank')
            end

            % Find Cz in channels
            elecLoc.(participantStr) = strcmpi(electrode,{EEG.chanlocs.labels});

            %% Extract relevant events of walking conditioons

            % Find if steps were taken during trial
            if sum(ismember([EEG.urevent(1:end).type], Event(Event == 100))) ~= 0

                % Find step onsets
                for iEvent = 1:2 %Correspond to uncued step onsets (preferred cadence)
                    j = 1;
                    for iLatency = 1:length(EEG.urevent)
                        if EEG.urevent(iLatency).type == Event(iEvent)
                            stepOnset(iEvent,j) = EEG.urevent(iLatency).latency;
                            j = j+1;
                        end
                    end
                end

                % Right foot
                stepOnset(stepOnset == 0) = NaN;
                stepR = diff(stepOnset(1,:));
                locsR = [];
                for iStepR = 1:length(stepR)
                    if stepR(iStepR) > 5000
                        locsR = [locsR; iStepR];
                    end
                end
                % figure; plot(stepR); hold on; plot(locsR,(stepR(locsR)), 'r*')
                eventLatency(1,:) = stepOnset(1,locsR); %Last right step onset before cueing

                
                eventLatency(1,end+1)  = max(stepOnset(1,:));          %Very last right step onset before cueing
                stepBegin(1,1) = stepOnset(1,1);                       %Very first uncued right step onset
                stepBegin(1,2:length(locsR)+1) = stepOnset(1,locsR+1); %First uncued right step onset after tempo shift

                % Left foot
                stepL = diff(stepOnset(2,:));
                locsL = [];
                for iStepL = 1:length(stepL)
                    if stepL(iStepL) > 5000
                        locsL = [locsL; iStepL];
                    end
                end
                % figure; plot(stepL); hold on; plot(locsL,(stepL(locsL)), 'r*')
                eventLatency(2,1:length(stepOnset(2,locsL))) = stepOnset(2,locsL); %Last left step onset before tempo cueing

                
                eventLatency(2,end+1)  = max(stepOnset(2,:));          %Very last left step onset before cueing
                stepBegin(2,1) = stepOnset(2,1);                       %Very first uncued right step onset
                stepBegin(2,2:length(locsL)+1) = stepOnset(2,locsL+1); %First uncued right step onset after tempo shift

                % Find important markers
                for iEvent = 3:length(Event)
                    j = 1;
                    for iLatency = 1:length(EEG.urevent)
                        if EEG.urevent(iLatency).type == Event(iEvent)
                            eventLatency(iEvent,j) = EEG.urevent(iLatency).latency;
                            j = j+1;
                        end
                    end
                end
                eventLatency(eventLatency == 0) = NaN;
                nCues = length(eventLatency(4,:)) + length(eventLatency(5,:));
                eventCues(1,:)  = sort(reshape(eventLatency(4:5,:),   [1, nCues]));
                eventCues(2,:)  = sort(reshape(eventLatency(11:12,:), [1, nCues]));
                eventCues(3,:)  = sort(reshape(eventLatency(13:14,:), [1, nCues]));
                eventCues(4,:)  = sort(reshape(eventLatency(1:2,:),   [1, nCues]));
                
                eventLatency      = eventLatency(3:6,:);
                stepBegin         = sort(reshape(stepBegin, [1 length(stepBegin)*2])); %Pool all right and left foot onsets of uncued walking together
                eventLatency(4,:) = nan(1,size(eventLatency,2));
                eventLatency(4,1:length(stepBegin(1:2:end))) = stepBegin(1:2:end); % Only keep one first heel strike (left of right) for uncued walking periods

                for iCond = 1:length(Conditions)-1 %Pref: iCond = 1; Slow: iCond = 2; Fast: iCond = 3; Uncued Pref: iCond = 4; cuedPrefRest: iCond = 5

                    if iBlock == 1
                        i = 1; %Total number of trials (all blocks)
                    else
                        i = size((dataAll.(participantStr).([Conditions{iCond}])),3) + 1;
                    end

                    for iTrial = 1:length(eventLatency(~isnan(eventLatency(iCond,:))))

                        if ~isnan(eventLatency(iCond,iTrial))

                            % Find time periods of cued gait at preferred cadence
                            if iCond == 1
                                if ~isnan(eventCues(iCond, iTrial))
                                    timeWin(iCond,:) = [eventLatency(iCond,iTrial) eventCues(iCond, iTrial)-1];
                                elseif isnan(eventCues(iCond, iTrial)) && ismember(EEG.urevent(end).type, [101 151 201 251])
                                    timeWin(iCond,:) = [eventLatency(iCond,iTrial) EEG.urevent(end).latency];
                                end

                            % Find time periods of cued gait at slow and fast cadences
                            else
                                target = eventCues(iCond,:) - eventLatency(iCond,iTrial+1);
                                target = target(target<0);
                                if ~isempty(target)
                                    timeWin(iCond,:) = [eventLatency(iCond,iTrial) eventCues(iCond, length(target))];
                                else
                                    Index = ~isnan(eventCues(iCond,:));
                                    timeWin(iCond,:) = [eventLatency(iCond,iTrial)  eventCues(iCond,sum(Index))];
                                end
                            end

                            if i == 1
                                dataAll.(participantStr).([Conditions{iCond}])(:,:,i) = nan(nChan,diff(timeWin(iCond,:))+1);
                            else
                                timeDiff = diff(timeWin(iCond,:))+1 - size((dataAll.(participantStr).([Conditions{iCond}])),2);

                                if timeDiff > 0
                                    dataAll.(participantStr).([Conditions{iCond}])(:,size((dataAll.(participantStr).([Conditions{iCond}])),2)+1:size((dataAll.(participantStr).([Conditions{iCond}])),2)+timeDiff,1:i-1) = nan;
                                    dataAll.(participantStr).([Conditions{iCond}])(:,:,i) = nan(nChan,diff(timeWin(iCond,:))+1);
                                elseif timeDiff <= 0
                                    dataAll.(participantStr).([Conditions{iCond}])(:,:,i) = nan(nChan,size((dataAll.(participantStr).([Conditions{iCond}])),2));
                                end

                            end
                            dataAll.(participantStr).([Conditions{iCond}])(:,1:diff(timeWin(iCond,:))+1,i) = data(:,timeWin(iCond,1):timeWin(iCond,2));
                            dataLengthTemp(i,iCond) = diff(timeWin(iCond,:))+1;

                            % Compute step freq
                            stepPrefCued   = [];
                            stepSlowCued   = [];
                            stepFastCued   = [];
                            stepPrefUncued = [];
                            for iLatency = 1:length(EEG.urevent)
                                if EEG.urevent(iLatency).latency >= timeWin(iCond,1) && EEG.urevent(iLatency).latency <= timeWin(iCond,2)
                                    if EEG.urevent(iLatency).type == 101     % Preferred cadence (right heel strike)
                                        stepPrefCued = [stepPrefCued; EEG.urevent(iLatency).latency];
                                    elseif EEG.urevent(iLatency).type == 201 % Preferred cadence (left heel strike)
                                        stepPrefCued = [stepPrefCued; EEG.urevent(iLatency).latency];
                                    elseif EEG.urevent(iLatency).type == 109 % Slow cadence (right heel strike)
                                        stepSlowCued = [stepSlowCued; EEG.urevent(iLatency).latency];
                                    elseif EEG.urevent(iLatency).type == 209 % Slow cadence (left heel strike)
                                        stepSlowCued = [stepSlowCued; EEG.urevent(iLatency).latency];
                                    elseif EEG.urevent(iLatency).type == 103 % Fast cadence (right heel strike)
                                        stepFastCued = [stepFastCued; EEG.urevent(iLatency).latency];
                                    elseif EEG.urevent(iLatency).type == 203 % Fast cadence (left heel strike)
                                        stepFastCued = [stepFastCued; EEG.urevent(iLatency).latency];
                                    elseif EEG.urevent(iLatency).type == 100 % Uncued preferred cadence (right heel strike)
                                        stepPrefUncued = [stepPrefUncued; EEG.urevent(iLatency).latency];
                                    elseif EEG.urevent(iLatency).type == 200 % Uncued preferred cadence (left heel strike)
                                        stepPrefUncued = [stepPrefUncued; EEG.urevent(iLatency).latency];                
                                    end
                                end
                            end

                            if iCond == 1
                                stepPrefCued = sort(stepPrefCued);
                                Freq(i,iCond) = mean(diff(stepPrefCued));
                            elseif iCond == 2
                                stepSlowCued = sort(stepSlowCued);
                                Freq(i,iCond) = mean(diff(stepSlowCued));
                            elseif iCond == 3
                                stepFastCued = sort(stepFastCued);
                                Freq(i,iCond) = mean(diff(stepFastCued));
                            elseif iCond == 4
                                stepPrefUncued = sort(stepPrefUncued);
                                Freq(i,iCond) = mean(diff(stepPrefUncued));
                            end
                            Freq(i,iCond) = rateEEG.(participantStr)/Freq(i,iCond);

                            i = i + 1;
                            clear stepPrefCued stepSlowCued stepFastCued target

                        end

                    end

                end

            % If there is no steps in trial (ie., standing condition)
            elseif sum(ismember([EEG.urevent(1:end).type], Event(Event == 171))) ~= 0
                iCond = 5;

                restPrefCued = [];
                if sum(ismember([EEG.urevent(1:end).type], [171 271])) == length(EEG.urevent)
                    for iLatency = 1:length(EEG.urevent)
                        restPrefCued = [restPrefCued; EEG.urevent(iLatency).latency];
                    end
                    
                    % Finding time periods of standing data
                    timeWin(iCond,1) = restPrefCued(1);
                    timeWin(iCond,2) = restPrefCued(end);

                    % Extracting rest data
                    if k == 1
                        dataAll.(participantStr).([Conditions{iCond}])(:,:,k) = nan(nChan,diff(timeWin(iCond,:))+1);
                    else
                        
                        timeDiff = diff(timeWin(iCond,:))+1 - size((dataAll.(participantStr).([Conditions{iCond}])),2);
                        if timeDiff > 0
                            dataAll.(participantStr).([Conditions{iCond}])(:,size((dataAll.(participantStr).([Conditions{iCond}])),2)+1:size((dataAll.(participantStr).([Conditions{iCond}])),2)+timeDiff,1:k-1) = nan;
                            dataAll.(participantStr).([Conditions{iCond}])(:,:,k) = nan(nChan,diff(timeWin(iCond,:))+1);
                        elseif timeDiff <= 0
                            dataAll.(participantStr).([Conditions{iCond}])(:,:,k) = nan(nChan,size((dataAll.(participantStr).([Conditions{iCond}])),2));
                        end

                    end
                    dataAll.(participantStr).([Conditions{iCond}])(:,1:diff(timeWin(iCond,:))+1,k) = data(:,timeWin(iCond,1):timeWin(iCond,2));
                    dataLengthTemp(k,iCond) = diff(timeWin(iCond,:))+1;

                    Freq(k,iCond) = mean(diff(restPrefCued));
                    Freq(k,iCond) = rateEEG.(participantStr)/Freq(k,iCond);
                    k = k + 1;

                    clear restPrefCued
                else
                    a = 1;
                end
                
            end

            ALLEEG = [];
            clear data eventLatency eventCues locsL locsR stepL stepR stepOnset stepRealR stepRealL stepBegin

        end

    end
     
    Freq(Freq == 0) = nan;
    freqAll.(participantStr) = Freq;

    dataLengthTemp(dataLengthTemp == 0) = nan;
    dataLength.(participantStr) = dataLengthTemp;

    save([pathResults 'subAll/DATA'], 'dataAll', 'dataLength', 'freqAll', 'elecLoc', 'rateEEG', 'chanLocs');

    clear dataLengthTemp Freq EEG 

end