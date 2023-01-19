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

% Load EEG variables computed with previous script
load([pathResults 'resultsEEG.mat'])

Participants = {'SUB001'; 'SUB002'; 'SUB003'; 'SUB004'; 'SUB005'; 'SUB006'; 'SUB007'; 'SUB008'};
Conditions   = {'Pref'; 'Slow'; 'Fast'; 'uncuedPref'; 'cuedPrefRest'};

% Pre-allocate matrices
power          = nan(length(Participants),length(Conditions));
stabilityIndex = nan(length(Participants),length(Conditions));

% Figure parameters
Colors     = [rgb('DodgerBlue'); rgb('DodgerBlue'); rgb('DodgerBlue'); rgb('DodgerBlue'); rgb('DodgerBlue')];
Titles     = {'Preferred Walk'; 'Slow Walk'; 'Fast Walk'; 'Preferred Walk (Control)'; 'Rest'};
titleFigs  = {'preferredWalk'; 'slowWalk'; 'fastWalk'; 'uncuedWalk'; 'preferredRest'};

eeglab;
for iConditions = 1:length(Conditions)

    for iParticipant = 1:length(Participants)
        
        %% Data extraction

        Freq     = resultsEEG.([Participants{iParticipant}]).([Conditions{iConditions}]).freqStep;
        EEGrate  = resultsEEG.([Participants{iParticipant}]).([Conditions{iConditions}]).freqEEG;
        ress     = resultsEEG.([Participants{iParticipant}]).([Conditions{iConditions}]).compTime';
        ressSNR  = resultsEEG.([Participants{iParticipant}]).([Conditions{iConditions}]).compSNR';
        chanLocs = resultsEEG.([Participants{iParticipant}]).([Conditions{iConditions}]).chanLocs;
        map      = resultsEEG.([Participants{iParticipant}]).([Conditions{iConditions}]).Map;

        %% Topography
        
        figure(iConditions+1);
        subplot(1,length(Participants), iParticipant);...
            topoplot(map./max(map), chanLocs, 'maplimits', [-1 1], 'numcontour', 0, 'conv', 'off', 'electrodes', 'on', 'shading', 'interp'); hold on;
            title(Participants{iParticipant})

        %% Power

        nFFT      = ceil(EEGrate/.02); % FFT resolution of .02 Hz
        Hz        = linspace(0,EEGrate,nFFT);
        freqIndex = dsearchn(Hz', Freq);        
        power(iParticipant, iConditions) = max(ressSNR(freqIndex-5:freqIndex+5));
     
        %% Stability index
        stabilityIndex(iParticipant,iConditions) = resultsEEG.([Participants{iParticipant}]).([Conditions{iConditions}]).stabilityIndex;

        clear ress chanLocs map
    end

end

%% Stats

% ANOVA
% tSI = table(Titles,stabilityIndex(:,1),stabilityIndex(:,2),stabilityIndex(:,3),stabilityIndex(:,4), stabilityIndex(:,5),...
% 'VariableNames',{'species','meas1','meas2','meas3','meas4'});
% Meas = table([1 2 3 4]','VariableNames',{'Measurements'});

% Power delta (SNR during cued gait minus SNR during uncued gait at preferred cadence)
powerDelta = power(:,1) - power(:,4);
[powerDP, powerDH, powerDStats] = signrank(powerDelta);

% Stability index delta (stability index during cued gait minus stability index during uncued gait at preferred cadence)
siDelta    = stabilityIndex(:,1) - stabilityIndex(:,4);
[siDP, siDH, siDStats]          = signrank(siDelta);

%% Plot

% Power
figure(7);
UnivarScatter(power, 'PointSize', 400, 'LineWidth', 1.5, 'MarkerFaceColor', Colors(:,:),'MarkerEdgeColor', 'none', 'Whiskers', 'lines', 'Compression', 2); hold on;
ax1 = gca;
xlim = get(ax1, 'xlim');
set(ax1, 'xticklabel', Titles);
ylabel('SNR')
set(ax1, 'FontWeight', 'bold', 'FontSize', 16);
title('Power')

% Power Delta
figure(8);
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

% Stability Index
figure(9);
UnivarScatter(stabilityIndex, 'PointSize', 400, 'LineWidth', 1.5, 'MarkerFaceColor', Colors(:,:),'MarkerEdgeColor', 'none', 'Whiskers', 'lines', 'Compression', 2); hold on;
ax3 = gca;
xlim = get(ax3, 'xlim');
set(ax3, 'xticklabel', Titles);
ylabel('Stability Index')
set(ax3, 'FontWeight', 'bold', 'FontSize', 16);
title('Stability Index')

% Stability Index Delta
figure(10);
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
saveas(figure(7),  [pathResults '/fig_Power.png']);
saveas(figure(8),  [pathResults '/fig_powerDelta.png']);
saveas(figure(9),  [pathResults '/fig_stabilityIndex.png']);
saveas(figure(10), [pathResults '/fig_stabilityIndexDelta.png']);