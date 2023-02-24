%% Compare power and stability index for all coniditions
% 1. Extract all relevant data from resultsEEG matrix
% 2. Insepct topographies for all participants and conditions
% 3. Compare power between all conditions
% 4. Compare stability index betweeen all conditions
% 5. Plot and save plots

clear all;
close all;
clc;

% Declare paths
pathResults = '/Volumes/Seagate/project_rhythmicBrain/Results/subAll/'; %Results path
addpath('/Users/claraziane/Documents/Académique/Informatique/MATLAB/eeglab2021.1'); %EEGLAB toolbox
addpath('/Volumes/Seagate/project_rhythmicBrain/Toolbox/rgb'); %To draw figures
addpath('/Volumes/Seagate/project_rhythmicBrain/Toolbox/UnivarScatter-master/UnivarScatter-master'); %To draw figures
addpath('/Users/claraziane/OneDrive - Universite de Montreal/S2M/projetRAC/dataAnalysis/Toolbox/CircStat2012a/'); %To compute resultant vector lengths

% Load EEG variables computed with previous script
load([pathResults 'resultsEEG.mat'])

Participants = {'SUB001'; 'SUB002'; 'SUB003'; 'SUB004'; 'SUB005'; 'SUB006'; 'SUB007'; 'SUB008'; 'SUB009'; 'SUB010'; 'SUB011'; 'SUB012'; 'SUB013'; 'SUB014'; 'SUB015'; 'SUB016'; 'SUB017'; 'SUB018'; 'SUB019'; 'SUB020'};
Conditions   = {'Pref'; 'Slow'; 'Fast'; 'uncuedPref'; 'cuedPrefRest'};
Variables    = {'stabilityIndex'; 'power'; 'stepR'; 'beatR'};

% Pre-allocate matrices
power          = nan(length(Participants),length(Conditions));
stabilityIndex = nan(length(Participants),length(Conditions));
stepR          = nan(length(Participants),length(Conditions)-1);
beatR          = nan(length(Participants),length(Conditions)-1);

% Figure parameters
Colors     = [rgb('DodgerBlue'); rgb('DodgerBlue'); rgb('DodgerBlue'); rgb('DodgerBlue'); rgb('DodgerBlue')];
Titles     = {'Preferred Walk'; 'Slow Walk'; 'Fast Walk'; 'Preferred Walk (Control)'; 'Rest'};
titleFigs  = {'preferredWalk'; 'slowWalk'; 'fastWalk'; 'uncuedWalk'; 'preferredRest'};
ylabels    = {'Stability Index'; 'SNR'; 'Resultant Vector Length'; 'Resultant Vector Length'};

eeglab;
for iCondition = 1:length(Conditions)

    for iParticipant = 1:length(Participants)
        
        %% Data extraction

        Freq      = resultsEEG.([Participants{iParticipant}]).([Conditions{iCondition}]).freqMax;
        EEGrate   = resultsEEG.([Participants{iParticipant}]).([Conditions{iCondition}]).freqEEG;
%         ress     = resultsEEG.([Participants{iParticipant}]).([Conditions{iConditions}]).compTime';
        ressSNR   = resultsEEG.([Participants{iParticipant}]).([Conditions{iCondition}]).compSNR';
        chanLocs  = resultsEEG.([Participants{iParticipant}]).([Conditions{iCondition}]).chanLocs;
        map       = resultsEEG.([Participants{iParticipant}]).([Conditions{iCondition}]).Map;
        stepPhase = resultsEEG.([Participants{iParticipant}]).([Conditions{iCondition}]).stepPhase;
        beatPhase = resultsEEG.([Participants{iParticipant}]).([Conditions{iCondition}]).beatPhase;

        %% Topography
        
        figure(iCondition+1);
        subplot(1,length(Participants), iParticipant);...
            topoplot(mean(map,2)./max(mean(map,2)), chanLocs, 'maplimits', [-1 1], 'numcontour', 0, 'conv', 'off', 'electrodes', 'on', 'shading', 'interp'); hold on;
            title(Participants{iParticipant})

        %% Power

        nFFT      = ceil(EEGrate/.02); % FFT resolution of .02 Hz
        Hz        = linspace(0,EEGrate,nFFT);
        for iBlock = 1:size(Freq,2)
            freqIndex = dsearchn(Hz', Freq(iBlock));   
            powerTemp(iBlock) = max(ressSNR(iBlock,freqIndex-5:freqIndex+5));
        end
        power(iParticipant, iCondition) = mean(powerTemp);
        resultsEEG.([Participants{iParticipant}]).([Conditions{iCondition}]).compPower = power(iParticipant, iCondition);
     
        %% Stability index
        stabilityIndex(iParticipant,iCondition) = resultsEEG.([Participants{iParticipant}]).([Conditions{iCondition}]).stabilityIndex;

        %% Phase
        if strcmpi(Conditions{iCondition}, 'cuedPrefRest') ~= 1
            stepR(iParticipant,iCondition) = circ_r(stepPhase, [], [], 1);
        end
        if strcmpi(Conditions{iCondition}, 'uncuedPref') ~= 1
            if strcmpi(Conditions{iCondition}, 'cuedPrefRest')
                beatR(iParticipant,iCondition-1) = circ_r(beatPhase, [], [], 1);
            else
                beatR(iParticipant,iCondition) = circ_r(beatPhase, [], [], 1);
            end
        end

        clear ress chanLocs map Freq powerTemp ressSNR stepPhase beatPhase
    end

end

% Save component power in results structure
save([pathResults 'resultsEEG'], 'resultsEEG');


%% Stats

% Power delta (SNR during cued gait minus SNR during uncued gait at preferred cadence)
powerDelta = power(:,1) - power(:,4);
[powerDP, powerDH, powerDStats] = signrank(powerDelta);

% Stability index delta (stability index during cued gait minus stability index during uncued gait at preferred cadence)
siDelta    = stabilityIndex(:,1) - stabilityIndex(:,4);
[siDP, siDH, siDStats]          = signrank(siDelta);

% Repeatd-measure ANOVA
for iVariable = 1:length(Variables) 

    if strcmpi(Variables{iVariable}, 'stepR')
        condTemp   = {'Pref'; 'Slow'; 'Fast'; 'uncuedPref'};
        titleTemp  = {'Preferred Walk'; 'Slow Walk'; 'Fast Walk'; 'Preferred Walk (Control)'};      
    elseif strcmpi(Variables{iVariable}, 'beatR')
        condTemp   = {'Pref'; 'Slow'; 'Fast'; 'cuedPrefRest'};
         titleTemp = {'Preferred Walk'; 'Slow Walk'; 'Fast Walk'; 'Rest'};      
    else
        condTemp = Conditions;
        titleTemp   = Titles;
    end

    % Create a data table
    data = array2table(eval(Variables{iVariable}), 'VariableNames', condTemp);

    % Indicate that "Conditions" is the within-subject factor
    factor = table(condTemp, 'VariableNames', {'Conditions'});

    % Fit the repeated-measure ANOVA model
    rm = fitrm(data, [condTemp{1} '-' condTemp{end} ' ~ 1'], 'WithinDesign', factor);

    % Run ANOVA
    ranovatbl = ranova(rm);

    % Tukey-Kramer post-hoc test
    Tukey = multcompare(rm, 'Conditions');

    %% Plot

    % Power
    figure(6+iVariable);
    UnivarScatter(eval(Variables{iVariable}), 'PointSize', 400, 'LineWidth', 1.5, 'MarkerFaceColor', Colors(:,:),'MarkerEdgeColor', 'none', 'Whiskers', 'lines', 'Compression', 2); hold on;
    ax = gca;
    xlim = get(ax, 'xlim');
    set(ax, 'xticklabel', titleTemp);
    ylabel(ylabels{iVariable})
    set(ax, 'FontWeight', 'bold', 'FontSize', 16);
    ylim = get(ax, 'ylim');
    i = 1;
    stop = 0;
    for iCompare = 1:length(Tukey.pValue)
        if Tukey.pValue(iCompare) < 0.05           
            tick1 = find(strcmpi(Tukey.Conditions_1(iCompare), condTemp) ~= 0);
            tick2 = find(strcmpi(Tukey.Conditions_2(iCompare), condTemp) ~= 0);   
            lineTicks(i,:) = sort(horzcat(tick1, tick2));
            for iLine = 1:size(lineTicks,1)
                if i > 1 && iLine < i ...
                    && ismember(lineTicks(i,1), lineTicks(iLine,1)) ...
                    && ismember(lineTicks(i,2), lineTicks(iLine,2))...
                    stop = 1;
                end
            end
            if stop == 1
                stop = 0;

            else
                set(ax, 'ylim', [ylim(1) ylim(2)+((ylim(2)-ylim(1))*0.10)]);
                ylim = get(ax, 'ylim');
                plot(lineTicks(i,:), [((ylim(2)-ylim(1))*0.95)+ylim(1) ((ylim(2)-ylim(1))*0.95)+ylim(1)], '-', 'color', 'black', 'LineWidth', 1.5,  'MarkerSize', 30);
                plot(min(horzcat(tick1,tick2))+abs(tick1 - tick2)/2, ((ylim(2)-ylim(1))*0.99)+ylim(1), '*', 'color', 'black', 'LineWidth', 1.5,  'MarkerSize', 15);
            end
            i = i + 1;
        end       
    end
    title(Variables{iVariable})
    clear lineTicks
    
    % Save figure
    saveas(figure(6+iVariable),  [pathResults '/fig_' Variables{iVariable} '.png']); 

end

% Power Delta
figure(11);
UnivarScatter(powerDelta, 'PointStyle', '^', 'PointSize', 400, 'LineWidth', 1.5, 'MarkerFaceColor', Colors(:,:),'MarkerEdgeColor', 'none', 'Whiskers', 'lines', 'Compression', 2); hold on;
ax2 = gca;
xlim = get(ax2, 'xlim');
plot([xlim(1) xlim(2)], [0 0], 'color', [0.80,0.80,0.80]);
set(ax2, 'xticklabel', 'Preferred Walk');
ylabel('SNR (\Delta)')
set(ax2, 'FontWeight', 'bold', 'FontSize', 16);
ylim = get(ax2, 'ylim');
if powerDH == 1
    set(ax2, 'ylim', [ylim(1) ylim(2)+((ylim(2)-ylim(1))*0.10)]);
    ylim = get(ax2, 'ylim');
    plot(((ylim(2)-ylim(1))*0.93)+ylim(1), '*', 'color', 'black', 'LineWidth', 1.5,  'MarkerSize', 30);
end
title('Power Delta')

% Stability Index Delta
figure(12);
UnivarScatter(siDelta, 'PointStyle', '^', 'PointSize', 400, 'LineWidth', 1.5, 'MarkerFaceColor', Colors(:,:),'MarkerEdgeColor', 'none', 'Whiskers', 'lines', 'Compression', 2); hold on;
ax4 = gca;
xlim = get(ax4, 'xlim');
plot([xlim(1) xlim(2)], [0 0], 'color', [0.80,0.80,0.80]);
set(ax4, 'xticklabel', 'Preferred Walk');
ylabel('Stability Index (\Delta)')
set(ax4, 'FontWeight', 'bold', 'FontSize', 16);
ylim = get(ax4, 'ylim');
if siDH == 1
    set(ax4, 'ylim', [ylim(1) ylim(2)+((ylim(2)-ylim(1))*0.10)]);
    ylim = get(ax4, 'ylim');
    plot(((ylim(2)-ylim(1))*0.93)+ylim(1), '*', 'color', 'black', 'LineWidth', 1.5,  'MarkerSize', 30);
end
title('Stability Index Delta')

% Save figures
saveas(figure(2),  [pathResults '/topo_' titleFigs{1} '.png']);
saveas(figure(3),  [pathResults '/topo_' titleFigs{2} '.png']);
saveas(figure(4),  [pathResults '/topo_' titleFigs{3} '.png']);
saveas(figure(5),  [pathResults '/topo_' titleFigs{4} '.png']);
saveas(figure(6),  [pathResults '/topo_' titleFigs{5} '.png']);
saveas(figure(11),  [pathResults '/fig_powerDelta.png']);
saveas(figure(12), [pathResults '/fig_stabilityIndexDelta.png']);

close all