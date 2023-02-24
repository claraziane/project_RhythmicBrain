clear all;
close all;
clc;

% Declare paths
pathData    = '/Volumes/Seagate/project_rhythmicBrain/DATA/';     %Folder where data was saved
pathResults = '/Volumes/Seagate/project_rhythmicBrain/Results/subAll/';  %Folder where to save results
addpath(genpath('/Volumes/Seagate/project_rhythmicBrain/Toolbox/rgb')); %To draw figures

Participants = {'sub-001'; 'sub-002'; 'sub-003'; 'sub-004'; 'sub-005'; 'sub-006'; 'sub-007'; 'sub-008'; 'sub-009'; 'sub-010'; 'sub-011'; 'sub-012'; 'sub-013'; 'sub-015'; 'sub-017'; 'sub-018'};
Conditions   = {'Pref'; 'Slow'; 'Fast'};
Colors       = [rgb('DodgerBlue'); rgb('DarkSlateBlue'); rgb('DarkOrange')];

% Load EEG & behavioural results previsouly computed
load([pathResults 'resultsBehavioural.mat']);
load([pathResults 'resultsEEG.mat']);

for iCondition = 1:length(Conditions)
        
    %% Extract relevant EEG and behavioural variables
    for iParticipant = 1:length(Participants)
        participantStr = strcat('SUB', Participants{iParticipant}(end-2:end));

        eegFreq = resultsEEG.(participantStr).([Conditions{iCondition}]).freqEEG;

        Power(iParticipant, iCondition)           = resultsEEG.(participantStr).([Conditions{iCondition}]).compPower;
        stabilityIndex(iParticipant, iCondition)  = resultsEEG.(participantStr).([Conditions{iCondition}]).stabilityIndex;

        IBIdeviation(iParticipant, iCondition)    = mean(resultsBehavioural.(participantStr).([Conditions{iCondition}]).IBIDeviation); %Compute block mean
        for iBlock = 1:size(eegFreq,1)
            asyncTemp(:,iBlock) = resultsBehavioural.(participantStr).([Conditions{iCondition}]).Asynchrony(:,iBlock) / (eegFreq/1000); %Convert to milliseconds
        end
        meanAsynchrony(iParticipant, iCondition)  = mean(mean(asyncTemp, 1));
        phaseAngleMean(iParticipant, iCondition)  = rad2deg(circ_mean(resultsBehavioural.(participantStr).([Conditions{iCondition}]).phaseAngleMean, [], 2));
        resultantLength(iParticipant,iCondition)  = nanmean(resultsBehavioural.(participantStr).([Conditions{iCondition}]).resultantLength);
        resultantLength(iParticipant,iCondition)  = log(resultantLength(iParticipant,iCondition)/(1-resultantLength(iParticipant,iCondition)));
        
        stepPhase = resultsEEG.(participantStr).([Conditions{iCondition}]).stepPhase;
        beatPhase = resultsEEG.(participantStr).([Conditions{iCondition}]).beatPhase;
        stepR(iParticipant,iCondition) = circ_r(stepPhase, [], [], 1);
        beatR(iParticipant,iCondition) = circ_r(beatPhase, [], [], 1);


        clear asyncTemp eegFreq stepPhase beatPhase

    end
    
    %% Compute correlations

    % Compute correlations between power and behavioural variables of synchronization
    [rhoIBI_P(iCondition),pIBI_P(iCondition)]     = corr(Power(:,iCondition), IBIdeviation(:, iCondition),    'Type', 'Spearman');
    [rhoAsync_P(iCondition),pAsync_P(iCondition)] = corr(Power(:,iCondition), meanAsynchrony(:,iCondition),   'Type', 'Spearman');
    [rhoPA_P(iCondition),pPA_P(iCondition)]       = corr(Power(:,iCondition), phaseAngleMean(:,iCondition),   'Type', 'Spearman');
    [rhoRVL_P(iCondition),pRVL_P(iCondition)]     = corr(Power(:,iCondition), resultantLength(:,iCondition),  'Type', 'Spearman');


    % Compute correlations between stability index and behavioural variables of synchronization
    [rhoIBI_SI(iCondition),pIBI_SI(iCondition)]     = corr(stabilityIndex(:,iCondition), IBIdeviation(:, iCondition),    'Type', 'Spearman');
    [rhoAsync_SI(iCondition),pAsync_SI(iCondition)] = corr(stabilityIndex(:,iCondition), meanAsynchrony(:,iCondition),   'Type', 'Spearman');
    [rhoPA_SI(iCondition),pPA_SI(iCondition)]       = corr(stabilityIndex(:,iCondition), phaseAngleMean(:,iCondition),   'Type', 'Spearman');
    [rhoRVL_SI(iCondition),pRVL_SI(iCondition)]     = corr(stabilityIndex(:,iCondition), resultantLength(:,iCondition),  'Type', 'Spearman');

    % Compute correlations between step RVL and behavioural variables of synchronization
    [rhoIBI_stepR(iCondition),pIBI_stepR(iCondition)]     = corr(stepR(:,iCondition), IBIdeviation(:, iCondition),    'Type', 'Spearman');
    [rhoAsync_stepR(iCondition),pAsync_stepR(iCondition)] = corr(stepR(:,iCondition), meanAsynchrony(:,iCondition),   'Type', 'Spearman');
    [rhoPA_stepR(iCondition),pPA_stepR(iCondition)]       = corr(stepR(:,iCondition), phaseAngleMean(:,iCondition),   'Type', 'Spearman');
    [rhoRVL_stepR(iCondition),pRVL_stepR(iCondition)]     = corr(stepR(:,iCondition), resultantLength(:,iCondition),  'Type', 'Spearman');

    % Compute correlations between beat RVL and behavioural variables of synchronization
    [rhoIBI_beatR(iCondition),pIBI_beatR(iCondition)]     = corr(beatR(:,iCondition), IBIdeviation(:, iCondition),    'Type', 'Spearman');
    [rhoAsync_beatR(iCondition),pAsync_beatR(iCondition)] = corr(beatR(:,iCondition), meanAsynchrony(:,iCondition),   'Type', 'Spearman');
    [rhoPA_beatR(iCondition),pPA_beatR(iCondition)]       = corr(beatR(:,iCondition), phaseAngleMean(:,iCondition),   'Type', 'Spearman');
    [rhoRVL_beatR(iCondition),pRVL_beatR(iCondition)]     = corr(beatR(:,iCondition), resultantLength(:,iCondition),  'Type', 'Spearman');

    % Pool all conditions together
    if iCondition == length(Conditions)
        PowerAll           = reshape(Power,           [length(Participants)*length(Conditions) 1]);
        stabilityIndexAll  = reshape(stabilityIndex,  [length(Participants)*length(Conditions) 1]);
        IBIdeviationAll    = reshape(IBIdeviation,    [length(Participants)*length(Conditions) 1]);
        meanAsynchronyAll  = reshape(meanAsynchrony,  [length(Participants)*length(Conditions) 1]);
        phaseAngleMeanAll  = reshape(phaseAngleMean,  [length(Participants)*length(Conditions) 1]);
        resultantLengthAll = reshape(resultantLength, [length(Participants)*length(Conditions) 1]);
        stepRAll           = reshape(stepR,           [length(Participants)*length(Conditions) 1]);
        beatRAll           = reshape(beatR,           [length(Participants)*length(Conditions) 1]);

        % Compute pooled correlations for power
        [rhoIBI_P_all,pIBI_P_all]     = corr(PowerAll, IBIdeviationAll,    'Type', 'Spearman');
        [rhoAsync_P_all,pAsync_P_all] = corr(PowerAll, meanAsynchronyAll,  'Type', 'Spearman');
        [rhoPA_P_all,pPA_P_all]       = corr(PowerAll, phaseAngleMeanAll,  'Type', 'Spearman');
        [rhoRVL_P_all,pRVL_P_all]     = corr(PowerAll, resultantLengthAll, 'Type', 'Spearman');

        % Compute pooled correlations for stability index
        [rhoIBI_SI_all,pIBI_SI_all]     = corr(stabilityIndexAll, IBIdeviationAll,    'Type', 'Spearman');
        [rhoAsync_SI_all,pAsync_SI_all] = corr(stabilityIndexAll, meanAsynchronyAll,  'Type', 'Spearman');
        [rhoPA_SI_all,pPA_SI_all]       = corr(stabilityIndexAll, phaseAngleMeanAll,  'Type', 'Spearman');
        [rhoRVL_SI_all,pRVL_SI_all]     = corr(stabilityIndexAll, resultantLengthAll, 'Type', 'Spearman');

        % Compute pooled correlations between step RVL and behavioural variables of synchronization
        [rhoIBI_stepR_all,pIBI_stepR_all]     = corr(stepRAll, IBIdeviationAll,    'Type', 'Spearman');
        [rhoAsync_stepR_all,pAsync_stepR_all] = corr(stepRAll, meanAsynchronyAll,   'Type', 'Spearman');
        [rhoPA_stepR_all,pPA_stepR_all]       = corr(stepRAll, phaseAngleMeanAll,   'Type', 'Spearman');
        [rhoRVL_stepR_all,pRVL_stepR_all]     = corr(stepRAll, resultantLengthAll,  'Type', 'Spearman');

        % Compute pooled correlations between beat RVL and behavioural variables of synchronization
        [rhoIBI_beatR_all,pIBI_beatR_all]     = corr(beatRAll, IBIdeviationAll,    'Type', 'Spearman');
        [rhoAsync_beatR_all,pAsync_beatR_all] = corr(beatRAll, meanAsynchronyAll,   'Type', 'Spearman');
        [rhoPA_beatR_all,pPA_beatR_all]       = corr(beatRAll, phaseAngleMeanAll,   'Type', 'Spearman');
        [rhoRVL_beatR_all,pRVL_beatR_all]     = corr(beatRAll, resultantLengthAll,  'Type', 'Spearman');

    end

    %% Plot

    % Plot correlations with power
    g = figure(1);
    sgtitle('Auditory-Motor Synchronization as a Function of Power', 'FontSize', 20)  
    ax = gca;    
    subplot(2,2,1); scatter(Power(:,iCondition), meanAsynchrony(:,iCondition), 150, Colors(iCondition,:), 'filled', 'MarkerFaceAlpha', 0.7,...
        'DisplayName', ['\rho = ' num2str(round(rhoAsync_P(iCondition),2)) '; p = ' num2str(round(pAsync_P(iCondition),2))]);
        xlabel('SNR', 'FontSize', 20); ylabel('Mean Asynchrony (ms)', 'FontSize', 20);
        if iCondition == length(Conditions)
            title(['Mean Asynchrony' newline '\rho = ' num2str(round(rhoAsync_P_all,4)) '; p = ' num2str(round(pAsync_P_all,6))], 'FontSize', 20);
            ax.FontSize = 16;
        end
        legend('FontSize', 16);
        hold on;
    subplot(2,2,2); scatter(Power(:,iCondition), IBIdeviation(:,iCondition), 150, Colors(iCondition,:), 'filled', 'MarkerFaceAlpha', 0.7,...
        'DisplayName', ['\rho = ' num2str(round(rhoIBI_P(iCondition),2)) '; p = ' num2str(round(pIBI_P(iCondition),2))]);
        xlabel('SNR', 'FontSize', 20); ylabel('IBI Deviation', 'FontSize', 20);
        if iCondition == length(Conditions)
            title(['Interbeat Interval Deviation' newline '\rho = ' num2str(round(rhoIBI_P_all,4)) '; p = ' num2str(round(pIBI_P_all,4))], 'FontSize', 20);
            ax.FontSize = 16;
        end
        legend('FontSize', 16);
        hold on;
    subplot(2,2,3); scatter(Power(:,iCondition), phaseAngleMean(:,iCondition), 150, Colors(iCondition,:), 'filled', 'MarkerFaceAlpha', 0.7,...
        'DisplayName', ['\rho = ' num2str(round(rhoPA_P(iCondition),2)) '; p = ' num2str(round(pPA_P(iCondition),2))]);
        xlabel('SNR', 'FontSize', 20); ylabel('Relative Phase Angle (°)', 'FontSize', 20);
        if iCondition == length(Conditions)
            title(['Relative Phase Angle' newline '\rho = ' num2str(round(rhoPA_P_all,4)) '; p = ' num2str(round(pPA_P_all,10))], 'FontSize', 20);
            ax.FontSize = 16;
        end
        legend('FontSize', 16);
        hold on;
    subplot(2,2,4); scatter(Power(:,iCondition), resultantLength(:,iCondition), 150, Colors(iCondition,:), 'filled', 'MarkerFaceAlpha', 0.7,...
        'DisplayName', ['\rho = ' num2str(round(rhoRVL_P(iCondition),2)) '; p = ' num2str(round(pRVL_P(iCondition),2))]);
        xlabel('SNR', 'FontSize', 20); ylabel('Resultant Vector Length (logit)', 'FontSize', 20);
        if iCondition == length(Conditions)
            title(['Resultant Vector Length' newline '\rho = ' num2str(round(rhoRVL_P_all,4)) '; p = ' num2str(round(pRVL_P_all,4))], 'FontSize', 20);
            ax.FontSize = 16;
        end
        legend('FontSize', 16);
        hold on;    
    
    % Plot correlations with stability index
    h = figure(2);
    sgtitle('Auditory-Motor Synchronization as a Function of the Stability Index', 'FontSize', 20)  
    bx = gca;    
    subplot(2,2,1); scatter(stabilityIndex(:,iCondition), meanAsynchrony(:,iCondition), 150, Colors(iCondition,:), 'filled', 'MarkerFaceAlpha', 0.7,...
        'DisplayName', ['\rho = ' num2str(round(rhoAsync_SI(iCondition),2)) '; p = ' num2str(round(pAsync_SI(iCondition),2))]);
        xlabel('Stability Index (Hz)', 'FontSize', 20); ylabel('Mean Asynchrony (ms)', 'FontSize', 20);
        if iCondition == length(Conditions)
            title(['Mean Asynchrony' newline '\rho = ' num2str(round(rhoAsync_SI_all,4)) '; p = ' num2str(round(pAsync_SI_all,4))], 'FontSize', 20);
            bx.FontSize = 16;
        end
        legend('FontSize', 16);
        hold on;
    subplot(2,2,2); scatter(stabilityIndex(:,iCondition), IBIdeviation(:,iCondition), 150, Colors(iCondition,:), 'filled', 'MarkerFaceAlpha', 0.7,...
        'DisplayName', ['\rho = ' num2str(round(rhoIBI_SI(iCondition),2)) '; p = ' num2str(round(pIBI_SI(iCondition),2))]);
        xlabel('Stability Index (Hz)', 'FontSize', 20); ylabel('IBI Deviation', 'FontSize', 20);
        if iCondition == length(Conditions)
            title(['Interbeat Interval Deviation' newline '\rho = ' num2str(round(rhoIBI_SI_all,4)) '; p = ' num2str(round(pIBI_SI_all,4))], 'FontSize', 20);
            bx.FontSize = 16;
        end
        legend('FontSize', 16);
        hold on;
    subplot(2,2,3); scatter(stabilityIndex(:,iCondition), phaseAngleMean(:,iCondition), 150, Colors(iCondition,:), 'filled', 'MarkerFaceAlpha', 0.7,...
        'DisplayName', ['\rho = ' num2str(round(rhoPA_SI(iCondition),2)) '; p = ' num2str(round(pPA_SI(iCondition),2))]);
        xlabel('Stability Index (Hz)', 'FontSize', 20); ylabel('Relative Phase Angle (°)', 'FontSize', 20);
        if iCondition == length(Conditions)
            title(['Relative Phase Angle' newline '\rho = ' num2str(round(rhoPA_SI_all,4)) '; p = ' num2str(round(pPA_SI_all,4))], 'FontSize', 20);
            bx.FontSize = 16;
        end
        legend('FontSize', 16);
        hold on;
    subplot(2,2,4); scatter(stabilityIndex(:,iCondition), resultantLength(:,iCondition), 150, Colors(iCondition,:), 'filled', 'MarkerFaceAlpha', 0.7,...
        'DisplayName', ['\rho = ' num2str(round(rhoRVL_SI(iCondition),2)) '; p = ' num2str(round(pRVL_SI(iCondition),2))]);
        xlabel('Stability Index (Hz)', 'FontSize', 20); ylabel('Resultant Vector Length (logit)', 'FontSize', 20);
        if iCondition == length(Conditions)
            title(['Resultant Vector Length' newline '\rho = ' num2str(round(rhoRVL_SI_all,4)) '; p = ' num2str(round(pRVL_SI_all,4))], 'FontSize', 20);
            bx.FontSize = 16;
        end
        legend('FontSize', 16);
        hold on;    

    % Plot correlations with step phase
    k = figure(3);
    sgtitle('Auditory-Motor Synchronization as a Function of the Phase at Step Onset', 'FontSize', 20)  
    cx = gca;    
    subplot(2,2,1); scatter(stepR(:,iCondition), meanAsynchrony(:,iCondition), 150, Colors(iCondition,:), 'filled', 'MarkerFaceAlpha', 0.7,...
        'DisplayName', ['\rho = ' num2str(round(rhoAsync_stepR(iCondition),2)) '; p = ' num2str(round(pAsync_stepR(iCondition),2))]);
        xlabel('Resultant Vector Length', 'FontSize', 20); ylabel('Mean Asynchrony (ms)', 'FontSize', 20);
        if iCondition == length(Conditions)
            title(['Mean Asynchrony' newline '\rho = ' num2str(round(rhoAsync_stepR_all,4)) '; p = ' num2str(round(pAsync_stepR_all,4))], 'FontSize', 20);
            cx.FontSize = 16;
        end
        legend('FontSize', 16);
        hold on;
    subplot(2,2,2); scatter(stepR(:,iCondition), IBIdeviation(:,iCondition), 150, Colors(iCondition,:), 'filled', 'MarkerFaceAlpha', 0.7,...
        'DisplayName', ['\rho = ' num2str(round(rhoIBI_stepR(iCondition),2)) '; p = ' num2str(round(pIBI_stepR(iCondition),2))]);
        xlabel('Resultant Vector Length', 'FontSize', 20); ylabel('IBI Deviation', 'FontSize', 20);
        if iCondition == length(Conditions)
            title(['Interbeat Interval Deviation' newline '\rho = ' num2str(round(rhoIBI_stepR_all,4)) '; p = ' num2str(round(pIBI_stepR_all,4))], 'FontSize', 20);
            cx.FontSize = 16;
        end
        legend('FontSize', 16);
        hold on;
    subplot(2,2,3); scatter(stepR(:,iCondition), phaseAngleMean(:,iCondition), 150, Colors(iCondition,:), 'filled', 'MarkerFaceAlpha', 0.7,...
        'DisplayName', ['\rho = ' num2str(round(rhoPA_stepR(iCondition),2)) '; p = ' num2str(round(pPA_stepR(iCondition),2))]);
        xlabel('Resultant Vector Length', 'FontSize', 20); ylabel('Relative Phase Angle (°)', 'FontSize', 20);
        if iCondition == length(Conditions)
            title(['Relative Phase Angle' newline '\rho = ' num2str(round(rhoPA_stepR_all,4)) '; p = ' num2str(round(pPA_stepR_all,4))], 'FontSize', 20);
            cx.FontSize = 16;
        end
        legend('FontSize', 16);
        hold on;
    subplot(2,2,4); scatter(stepR(:,iCondition), resultantLength(:,iCondition), 150, Colors(iCondition,:), 'filled', 'MarkerFaceAlpha', 0.7,...
        'DisplayName', ['\rho = ' num2str(round(rhoRVL_stepR(iCondition),2)) '; p = ' num2str(round(pRVL_stepR(iCondition),2))]);
        xlabel('Resultant Vector Length', 'FontSize', 20); ylabel('Resultant Vector Length (logit)', 'FontSize', 20);
        if iCondition == length(Conditions)
            title(['Resultant Vector Length' newline '\rho = ' num2str(round(rhoRVL_stepR_all,4)) '; p = ' num2str(round(pRVL_stepR_all,4))], 'FontSize', 20);
            cx.FontSize = 16;
        end
        legend('FontSize', 16);
        hold on;    
        
    % Plot correlations with beat phase
    l = figure(4);
    sgtitle('Auditory-Motor Synchronization as a Function of the Phase at Beat Onset', 'FontSize', 20)  
    dx = gca;    
    subplot(2,2,1); scatter(beatR(:,iCondition), meanAsynchrony(:,iCondition), 150, Colors(iCondition,:), 'filled', 'MarkerFaceAlpha', 0.7,...
        'DisplayName', ['\rho = ' num2str(round(rhoAsync_beatR(iCondition),2)) '; p = ' num2str(round(pAsync_beatR(iCondition),2))]);
        xlabel('Resultant Vector Length', 'FontSize', 20); ylabel('Mean Asynchrony (ms)', 'FontSize', 20);
        if iCondition == length(Conditions)
            title(['Mean Asynchrony' newline '\rho = ' num2str(round(rhoAsync_beatR_all,4)) '; p = ' num2str(round(pAsync_beatR_all,4))], 'FontSize', 20);
            dx.FontSize = 16;
        end
        legend('FontSize', 16);
        hold on;
    subplot(2,2,2); scatter(beatR(:,iCondition), IBIdeviation(:,iCondition), 150, Colors(iCondition,:), 'filled', 'MarkerFaceAlpha', 0.7,...
        'DisplayName', ['\rho = ' num2str(round(rhoIBI_beatR(iCondition),2)) '; p = ' num2str(round(pIBI_beatR(iCondition),2))]);
        xlabel('Resultant Vector Length', 'FontSize', 20); ylabel('IBI Deviation', 'FontSize', 20);
        if iCondition == length(Conditions)
            title(['Interbeat Interval Deviation' newline '\rho = ' num2str(round(rhoIBI_beatR_all,4)) '; p = ' num2str(round(pIBI_beatR_all,4))], 'FontSize', 20);
            dx.FontSize = 16;
        end
        legend('FontSize', 16);
        hold on;
    subplot(2,2,3); scatter(beatR(:,iCondition), phaseAngleMean(:,iCondition), 150, Colors(iCondition,:), 'filled', 'MarkerFaceAlpha', 0.7,...
        'DisplayName', ['\rho = ' num2str(round(rhoPA_beatR(iCondition),2)) '; p = ' num2str(round(pPA_beatR(iCondition),2))]);
        xlabel('Resultant Vector Length', 'FontSize', 20); ylabel('Relative Phase Angle (°)', 'FontSize', 20);
        if iCondition == length(Conditions)
            title(['Relative Phase Angle' newline '\rho = ' num2str(round(rhoPA_beatR_all,4)) '; p = ' num2str(round(pPA_beatR_all,4))], 'FontSize', 20);
            dx.FontSize = 16;
        end
        legend('FontSize', 16);
        hold on;
    subplot(2,2,4); scatter(beatR(:,iCondition), resultantLength(:,iCondition), 150, Colors(iCondition,:), 'filled', 'MarkerFaceAlpha', 0.7,...
        'DisplayName', ['\rho = ' num2str(round(rhoRVL_beatR(iCondition),2)) '; p = ' num2str(round(pRVL_beatR(iCondition),2))]);
        xlabel('Resultant Vector Length', 'FontSize', 20); ylabel('Resultant Vector Length (logit)', 'FontSize', 20);
        if iCondition == length(Conditions)
            title(['Resultant Vector Length' newline '\rho = ' num2str(round(rhoRVL_beatR_all,4)) '; p = ' num2str(round(pRVL_beatR_all,4))], 'FontSize', 20);
            dx.FontSize = 16;
        end
        legend('FontSize', 16);
        hold on;    

end
saveas(figure(1), [pathResults 'fig_corrPower.png']);
saveas(figure(2), [pathResults 'fig_CorrSI.png']);
saveas(figure(3), [pathResults 'fig_CorrStepR.png']);
saveas(figure(4), [pathResults 'fig_CorrBeatR.png']);