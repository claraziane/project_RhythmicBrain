clear all;
close all;
clc;

% Declare paths
pathData    = '/Volumes/Seagate/project_rhythmicBrain/DATA/';     %Folder where data was saved
pathResults = '/Volumes/Seagate/project_rhythmicBrain/Results/subAll/';  %Folder where to save results
addpath('/Users/claraziane/OneDrive - Universite de Montreal/S2M/projetRAC/dataAnalysis/Toolbox/CircStat2012a/');

Participants = {'sub-001'; 'sub-002'; 'sub-003'; 'sub-004'; 'sub-005'; 'sub-006'; 'sub-007'; 'sub-008'; 'sub-009'; 'sub-010'; 'sub-011'; 'sub-012'; 'sub-013'; 'sub-015'; 'sub-017'; 'sub-018'};
Conditions   = {'Pref'; 'Slow'; 'Fast'};

% Load results already computed
% load([pathResults 'resultsBehavioural.mat']); %Comment if running code for the 1st time
load([pathResults 'resultsEEG.mat']); %Comment if running script for the 1st time

for iParticipant = 1:length(Participants)
    participantStr = strcat('SUB', Participants{iParticipant}(end-2:end));

    % Load data previously extracted for participant
    load([pathData Participants{iParticipant} '/eeg/equalDATA.mat'], 'cueData', 'stepData');

    for iCondition = 1:length(Conditions)

        % Extract trials to remove from EEG results
        trials2remove = resultsEEG.(participantStr).([Conditions{iCondition}]).trialsRemoved;
        iTrial = 1;

        % Pre-allocating matrices
        phaseAngle = nan(size(stepData.([Conditions{iCondition}]),1)-1, size(stepData.([Conditions{iCondition}]),2)-length(trials2remove));
        phaseRad   = nan(size(stepData.([Conditions{iCondition}]),1)-1, size(stepData.([Conditions{iCondition}]),2)-length(trials2remove));
        for iBlock = 1:size(stepData.([Conditions{iCondition}]),2)
            
            if ~ismember(iBlock, trials2remove)

                %% Estimating period-matching accuracy (i.e., extent to which step tempo matches stimulus tempo) using IBI deviation

                % Extracting beat onsets
                beatOnset = [];
                beatOnset = cueData.([Conditions{iCondition}])(:,iBlock);
                beatOnset = beatOnset(~isnan(beatOnset));

                % Extracting step onsets
                stepOnset = [];
                stepOnset = stepData.([Conditions{iCondition}])(:,iBlock);

                % Matching step onsets to closest beat
                beatMatched = [];
                for iStep = 1:length(stepOnset)
                    [minValue matchIndex] = min(abs(beatOnset-stepOnset(iStep)));
                    beatMatched(iStep,1) = beatOnset(matchIndex);
                end

                % Calculating interstep interval
                stepInterval = [];
                stepInterval = diff(stepOnset);

                % Calculating interbeat interval
                racInterval = [];
                racInterval = diff(beatMatched);

                % Calculating IBI deviation
                IBI(iTrial) = mean(abs(stepInterval - racInterval))/mean(racInterval);

                %% Estimating phase-matching accuracy (i.e., the difference between step onset times and beat onset times) using circular asynchronies

                asynchrony(:,iTrial)           = stepOnset - beatMatched;
                asynchronyNormalized(:,iTrial) = asynchrony(1:end-1,iTrial)./stepInterval;
                asynchronyCircular(:,iTrial)   = asynchronyNormalized(:,iTrial) * 360;
                asynchronyRad(:,iTrial)        = asynchronyCircular(:,iTrial) * pi/180;
                asynchronyMean(iTrial)         = circ_mean(asynchronyRad(:,iTrial), [], 1);
                %             figure; scatter(1,asynchronyCircular(:,iBlock))

                % Running Rao's test (a not-significant test means participant failed to synchronize)
                [p(iTrial) U UC] = circ_raotest(asynchronyCircular(:,iTrial));

                % Calculating circular variance
                [varianceCircular(iTrial) varianceAngular(iTrial)] = circ_var(asynchronyCircular(:,iTrial)/(180/pi)); % Degrees are converted to radians

                % Calculating phase angles (error measure of synchronization based on the phase difference between two oscillators)
                phaseAngle(1:length(racInterval(racInterval ~= 0)),iTrial) = 360*((stepOnset(racInterval ~= 0) - beatMatched(racInterval ~= 0))./racInterval(racInterval ~= 0));
                phaseRad(1:length(racInterval(racInterval ~= 0)),iTrial)   = deg2rad(phaseAngle(~isnan(phaseAngle(:,iTrial)),iTrial));
                phaseAngleMean(iTrial) = circ_mean(phaseRad(~isnan(phaseRad(:,iTrial)),iTrial), [], 1);

                % Calculating resultant vector length (expresses the stability of the relative phase angles over time)
                resultantLength(iTrial) = circ_r(phaseRad(:,iTrial), [], [], 1);

                iTrial = iTrial +1;
            end

        end

        % Storing results in structure
        resultsBehavioural.(participantStr).([Conditions{iCondition}]).IBIDeviation = IBI;
        resultsBehavioural.(participantStr).([Conditions{iCondition}]).Asynchrony = asynchrony;
        resultsBehavioural.(participantStr).([Conditions{iCondition}]).circularAsynchrony = asynchronyCircular;
        resultsBehavioural.(participantStr).([Conditions{iCondition}]).asynchronyMean = asynchronyMean;
        resultsBehavioural.(participantStr).([Conditions{iCondition}]).circularVariance = varianceCircular;
        resultsBehavioural.(participantStr).([Conditions{iCondition}]).pRao = p;
        resultsBehavioural.(participantStr).([Conditions{iCondition}]).phaseAngle = phaseAngle;
        resultsBehavioural.(participantStr).([Conditions{iCondition}]).phaseAngleMean = phaseAngleMean;
        resultsBehavioural.(participantStr).([Conditions{iCondition}]).resultantLength = resultantLength;
        clear asynchrony asynchronyCircular asynchronyMean asynchronyNormalized  asynchronyRad ...
            IBI p phaseAngle phaseAngleMean phaseRad resultantLength VarianceAngular varianceCircular

    end

end

% Save results
save([pathResults 'resultsBehavioural.mat'], 'resultsBehavioural');