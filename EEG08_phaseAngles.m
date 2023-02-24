clear all;
close all;
clc;

% Declare paths
pathResults = '/Volumes/Seagate/project_rhythmicBrain/Results/'; %Results path
addpath('/Users/claraziane/OneDrive - Universite de Montreal/S2M/projetRAC/dataAnalysis/Toolbox/CircStat2012a/'); %To compute resultant vector lengths
addpath('/Users/claraziane/OneDrive - Universite de Montreal/S2M/projetRAC/dataAnalysis/Functions/circularHist');


% Load EEG variables computed with previous script
load([pathResults 'subAll/resultsEEG.mat'])

Participants = {'sub-001'; 'sub-002'; 'sub-003'; 'sub-004'; 'sub-005'; 'sub-006'; 'sub-007'; 'sub-008'; 'sub-009'; 'sub-010'; 'sub-011'; 'sub-012'; 'sub-013'; 'sub-015'; 'sub-017'; 'sub-018'};
Conditions   = {'Pref'; 'Slow'; 'Fast'; 'uncuedPref'; 'cuedPrefRest'};

for iParticipant = 1:length(Participants)
    participantStr = strcat('SUB', Participants{iParticipant}(end-2:end));

    for iCondition = 1:length(Conditions)

        if strcmpi(Conditions{iCondition}, 'cuedPrefRest') == 0
            phaseStep = resultsEEG.(participantStr).([Conditions{iCondition}]).stepPhase;
            phaseStep = phaseStep * (180/pi);

            % Draw phase angle figure
            nBins = 360/5;
            fH = figure(1);
            ax = polaraxes(fH);
            obj1 = CircHist(phaseStep, nBins, 'parent', ax);
            obj1.polarAxs.ThetaZeroLocation = 'right'; % rotate the plot to have 0° on the right side
            fH.Visible = 'on';
            obj1.polarAxs.Title.Parent.ThetaTickLabel = {'30°'; '60°'; '90°'; '120°'; '150°'; '180°'; '-150°'; '-120°'; '-90°'; '-60°'; '-30'; '0°'};
            obj1.thetaLabel.Parent.ThetaLim = [-359 0];

            % Save figures
            saveas(figure(1),  [pathResults Participants{iParticipant} '/Stable/' Conditions{iCondition} '/fig_phaseStep.png']);

        end

        if strcmpi(Conditions{iCondition}, 'uncuedPref') == 0
            phaseBeat = resultsEEG.(participantStr).([Conditions{iCondition}]).beatPhase;
            phaseBeat = phaseBeat * (180/pi);

            gH = figure(2);
            ax = polaraxes(gH);
            obj2 = CircHist(phaseBeat, nBins, 'parent', ax);
            obj2.polarAxs.ThetaZeroLocation = 'right'; % rotate the plot to have 0° on the right side
            gH.Visible = 'on';
            obj2.polarAxs.Title.Parent.ThetaTickLabel = {'30°'; '60°'; '90°'; '120°'; '150°'; '180°'; '-150°'; '-120°'; '-90°'; '-60°'; '-30'; '0°'};
            obj2.thetaLabel.Parent.ThetaLim = [-359 0];

            % Save figures
            saveas(figure(2),  [pathResults Participants{iParticipant} '/Stable/' Conditions{iCondition} '/fig_phaseBeat.png']);

        end

        close all;
        clear phaseStep phaseBeat

    end

end