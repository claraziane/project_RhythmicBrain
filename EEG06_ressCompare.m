clear all;
close all;
clc;

pathResults = '/Volumes/Seagate/project_rhythmicBrain/Results/subAll/';
addpath('/Users/claraziane/Documents/Académique/Informatique/MATLAB/eeglab2021.1');
% addpath '/Volumes/10.89.24.15/Projet_RAC/dataAnalysis/Scripts/Clara/EEG/Functions'
addpath('/Volumes/Seagate/project_rhythmicBrain/Toolbox/rgb'); %To draw figures
addpath('/Volumes/Seagate/project_rhythmicBrain/Toolbox/UnivarScatter-master/UnivarScatter-master'); %To draw figures
load('/Volumes/Seagate/project_rhythmicBrain/Results/subAll/resultsEEG.mat')
% load('/Volumes/10.89.24.15/Projet_RAC/DATA/Preprocessed/RAC/All/RAC.mat');
% load('/Volumes/10.89.24.15/Projet_RAC/DATA/Preprocessed/Kinetics/All/stepData.mat');
load([pathResults 'resultsEEG.mat'])

Participants = {'SUB001'; 'SUB002'; 'SUB003'; 'SUB004';};
Blocks       = {'Pref'; 'Slow'; 'Fast'; 'uncuedPref'; 'cuedPrefRest'};
power                 = nan(length(Participants),length(Blocks));
powerControl          = nan(length(Participants),length(Blocks));
stabilityIndex        = nan(length(Participants),length(Blocks));
stabilityIndexControl = nan(length(Participants),length(Blocks));

% Figure parameters
Colors     = [rgb('DodgerBlue'); rgb('DodgerBlue'); rgb('DodgerBlue'); rgb('DodgerBlue'); rgb('DodgerBlue')];
Titles     = {'Preferred Walk'; 'Slow Walk'; 'Fast Walk'; 'Preferred Walk (Control)'; 'Rest'};
titleFigs  = {'preferredWalk'; 'slowWalk'; 'fastWalk'; 'uncuedWalk'; 'preferredRest'};
% xLabels    = {'RAC'; 'Control'};

eeglab;
for iCompare = 1:length(Blocks)

    for iParticipant = 1:length(Participants)
        
        %% Data extraction
%         if strcmp(Blocks{1,iCompare}(end-3:end), 'Walk')
            Freq = resultsEEG.([Participants{iParticipant}]).([Blocks{iCompare}]).freqStep;
%             FreqControl = round(stepData.([Participants{iParticipant}]).([Blocks{2,iCompare}]).stepFreq,2);
%         elseif strcmp(Blocks{1,iCompare}(end-3:end), 'Rest')
%             Freq        = round(RAC.([Participants{iParticipant}]).([Blocks{1,iCompare}]).beatFrequency,2);
%             FreqControl = round(RAC.([Participants{iParticipant}]).([Blocks{1,iCompare}]).beatFrequency,2);
%         end

        EEGrate = resultsEEG.([Participants{iParticipant}]).([Blocks{iCompare}]).freqEEG;
        ress            = resultsEEG.([Participants{iParticipant}]).([Blocks{iCompare}]).compTime';
%         ressControl     = resultsEEG.([Participants{iParticipant}]).([Blocks{2,iCompare}]).compTime';
        ressSNR         = resultsEEG.([Participants{iParticipant}]).([Blocks{iCompare}]).compSNR';
%         controlSNR      = resultsEEG.([Participants{iParticipant}]).([Blocks{2,iCompare}]).compSNR';      
        chanLocs        = resultsEEG.([Participants{iParticipant}]).([Blocks{iCompare}]).chanLocs;
%         chanLocsControl = resultsEEG.([Participants{iParticipant}]).([Blocks{2,iCompare}]).chanLocs;
        map        = resultsEEG.([Participants{iParticipant}]).([Blocks{iCompare}]).Map;
%         mapControl = resultsEEG.([Participants{iParticipant}]).([Blocks{2,iCompare}]).Map;

        %% Topography
        figure(iCompare+1);
        subplot(1,length(Participants), iParticipant);...
            topoplot(map./max(map), chanLocs, 'maplimits', [-1 1], 'numcontour', 0, 'conv', 'off', 'electrodes', 'on', 'shading', 'interp'); hold on;
            title(Participants{iParticipant})
%         subplot(2, length(Participants), iParticipant+length(Participants));...
%             topoplot(mapControl, chanLocsControl, 'maplimits', [-0.7 0.7], 'numcontour', 0, 'conv', 'off', 'electrodes', 'on', 'shading', 'interp');

        %% Power

        nFFT = ceil(EEGrate/0.01);
        Hz    = linspace(0,EEGrate,nFFT);
        freqIndex        = dsearchn(Hz', Freq);
%         freqIndexControl = dsearchn(Hz', FreqControl);
        
        % RESS method
        power(iParticipant, iCompare)        = max(ressSNR(freqIndex-5:freqIndex+5));
%         powerControl(iParticipant, iCompare) = max(controlSNR(freqIndexControl-5:freqIndexControl+5));
     
        %% Stability index
        stabilityIndex(iParticipant,iCompare)        = resultsEEG.([Participants{iParticipant}]).([Blocks{iCompare}]).stabilityIndex;
%         stabilityIndexControl(iParticipant,iCompare) = resultsEEG.([Participants{iParticipant}]).([Blocks{2,iCompare}]).stabilityIndex;

        clear ress ressx ressControl controlx chanLocs map mapControl mapDiff
    end

    %% Stats
%     powerDelta(:,iCompare) = power(:,iCompare) - powerControl(:,iCompare);
%     siDelta(:,iCompare)    = stabilityIndex(:,iCompare) - stabilityIndexControl(:,iCompare);

%     [powerH(iCompare),  powerP(iCompare),  powerCI(iCompare,:),   powerStats] = ranksum(power(:,iCompare), powerControl(:, iCompare));
%     [ powerDP(iCompare), powerDH(iCompare),  powerDStats] = signrank(powerDelta(~isnan(powerDelta(:,iCompare)), iCompare)');
%     [siH(iCompare),     siP(iCompare),     siCI(iCompare,:),         siStats] = ranksum(stabilityIndex(:, iCompare), stabilityIndexControl(:, iCompare));
%     [siDP(iCompare), siDH(iCompare),             siDStats] = signrank(siDelta(~isnan(siDelta(:,iCompare)), iCompare));
    
    %% Plotting

%      % Power
%      figure(6);
%      subplot(1,length(Blocks),iCompare); fig6 = UnivarScatter([power(:,iCompare), powerControl(:,iCompare)], 'PointSize', 200, 'LineWidth', 1.5, 'MarkerFaceColor', Color(iCompare,:),'MarkerEdgeColor', 'none', 'Whiskers', 'lines', 'Compression', 2); hold on;
%      plot([fig6(:,1) fig6(:,2)]', [power(:,iCompare), powerControl(:,iCompare)]', 'Color', Colors(iCompare,:), 'LineWidth', 1, 'LineStyle', '-'); hold on;
%      ax1 = gca;
%      set(ax1, 'xticklabel', xLabels);
%      set(ax1, 'FontWeight', 'bold', 'FontSize', 16);
% %      if powerH(iCompare) == 1
% %          ylim = get(ax1, 'ylim');
% %          set(ax1, 'ylim', [ylim(1) ylim(2)+((ylim(2)-ylim(1))*0.10)]);
% %          ylim = get(ax1, 'ylim');
% %          plot(1.5, ((ylim(2)-ylim(1))*0.93)+ylim(1), '*', 'color', 'black', 'LineWidth', 1.5,  'MarkerSize', 30);
% %      end
%      if iCompare == 1
%          ylabel('SNR (Delta)')
%      end
%      title(Titles{iCompare})
% 
%      % Stability index
%      figure(7);
%      subplot(1,length(Blocks),iCompare); fig7 = UnivarScatter([stabilityIndex(:,iCompare), stabilityIndexControl(:,iCompare)], 'PointSize', 200, 'LineWidth', 1.5, 'MarkerFaceColor', Color(iCompare,:),'MarkerEdgeColor', 'none', 'Whiskers', 'lines','Compression', 2); hold on;
%      plot([fig7(:,1) fig7(:,2)]', [stabilityIndex(:,iCompare), stabilityIndexControl(:,iCompare)]', 'Color', Colors(iCompare,:), 'LineWidth', 1.5, 'LineStyle', '-'); hold on;
%      ax2 = gca;
%      set(ax2, 'xticklabel', xLabels);
%      set(ax2, 'FontWeight', 'bold', 'FontSize', 16);
% %      if siH(iCompare) == 1
% %          ylim = get(ax2, 'ylim');
% %          set(ax2, 'ylim', [ylim(1) ylim(2)+((ylim(2)-ylim(1))*0.10)]);
% %          ylim = get(ax2, 'ylim');
% %          plot(1.5, ((ylim(2)-ylim(1))*0.93)+ylim(1), '*', 'color', 'black', 'LineWidth', 1.5,  'MarkerSize', 30);
% %      end
%      if iCompare == 1
%          ylabel('Index de Stabilité')
%      end
%      title(Titles{iCompare})
end

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

% Stability Index
figure(8);
UnivarScatter(stabilityIndex, 'PointSize', 400, 'LineWidth', 1.5, 'MarkerFaceColor', Colors(:,:),'MarkerEdgeColor', 'none', 'Whiskers', 'lines', 'Compression', 2); hold on;
ax1 = gca;
xlim = get(ax1, 'xlim');
set(ax1, 'xticklabel', Titles);
ylabel('Stability Index')
set(ax1, 'FontWeight', 'bold', 'FontSize', 16);
title('Stability Index')

% Save figures
saveas(figure(2),  [pathResults '/topo_' titleFigs{1} '.png']);
saveas(figure(3),  [pathResults '/topo_' titleFigs{2} '.png']);
saveas(figure(4),  [pathResults '/topo_' titleFigs{3} '.png']);
saveas(figure(5),  [pathResults '/topo_' titleFigs{4} '.png']);
saveas(figure(6),  [pathResults '/topo_' titleFigs{5} '.png']);

% saveas(figure(5),  [pathResults '/fig_Power.png']);
% xsaveas(figure(6), [pathResults '/fig_stabilityIndex.png']);

% % Power
% figure(iCompare+3);
% UnivarScatter(powerDelta, 'PointStyle', '^', 'PointSize', 400, 'LineWidth', 1.5, 'MarkerFaceColor', Colors(:,:),'MarkerEdgeColor', 'none', 'Whiskers', 'lines', 'Compression', 2); hold on;
% ax3 = gca;
% xlim = get(ax3, 'xlim');
% plot([xlim(1) xlim(2)], [0 0], 'color', [0.80,0.80,0.80]);
% set(ax3, 'xticklabel', Titles);
% ylabel('SNR (\Delta)')
% set(ax3, 'FontWeight', 'bold', 'FontSize', 16);
% ylim = get(ax3, 'ylim');
% set(ax3, 'ylim', [ylim(1) ylim(2)+((ylim(2)-ylim(1))*0.10)]);
% ylim = get(ax3, 'ylim');
% for iCompare = 1:4
%     if powerDH(iCompare) == 1
%         plot(iCompare, ((ylim(2)-ylim(1))*0.93)+ylim(1), '*', 'color', 'black', 'LineWidth', 1.5,  'MarkerSize', 30);
%     end
% end
% title('Power')
% saveas(figure(iCompare+3), [pathResults analysisType '/fig_powerDelta.png']);
% 
% figure(iCompare+4);
% UnivarScatter(siDelta, 'PointStyle', '^', 'PointSize', 400, 'LineWidth', 1.5, 'MarkerFaceColor', Colors(:,:),'MarkerEdgeColor', 'none', 'Whiskers', 'lines', 'Compression', 2); hold on;
% ax4 = gca;
% xlim = get(ax4, 'xlim');
% plot([xlim(1) xlim(2)], [0 0], 'color', [0.80,0.80,0.80]);
% set(ax4, 'xticklabel', Titles);
% ylabel('Stability Index (\Delta)')
% set(ax4, 'FontWeight', 'bold', 'FontSize', 16);
% ylim = get(ax4, 'ylim');
% set(ax4, 'ylim', [ylim(1) ylim(2)+((ylim(2)-ylim(1))*0.10)]);
% ylim = get(ax4, 'ylim');
% for iCompare = 1:4
%     if siDH(iCompare) == 1
%         plot(iCompare, ((ylim(2)-ylim(1))*0.93)+ylim(1), '*', 'color', 'black', 'LineWidth', 1.5,  'MarkerSize', 30);
%     end
% end
% title('Stability Index')
% saveas(figure(iCompare+4), [pathResults analysisType '/fig_stabilityIndexDelta.png']);

%% Stats

tSI = table(Titles,stabilityIndex(:,1),stabilityIndex(:,2),stabilityIndex(:,3),stabilityIndex(:,4), stabilityIndex(:,5),...
'VariableNames',{'species','meas1','meas2','meas3','meas4'});
Meas = table([1 2 3 4]','VariableNames',{'Measurements'});