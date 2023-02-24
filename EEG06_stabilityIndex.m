%% Extracting entrained component & computing the stability index
% 1. Compute RESS component
%   a. sequence EEG in trials according to beat onsets (-100:500 ms)
%   b. compute covariance matrices S (from narrow-filtered data) and R (from neighbouring frequencies)
%   c. perform Generalized Eigen Decomposition (GED)
% 2. Compute stability index
%   a. Filter RESS component
%   b. Transform filtered RESS into analytical signal using Hilbert transform
%   c. Compute phase angles
%   d. Extract instantaneous frequencies
%   e. Compute stability index (standard deviation of instantaneous frequencies)

clear;
close all;
clc;

% Declare paths
pathData    = '/Volumes/Seagate/project_rhythmicBrain/DATA/';
pathResults = '/Volumes/Seagate/project_rhythmicBrain/Results/';  %Folder where to save results
addpath('/Users/claraziane/Documents/Académique/Informatique/MATLAB/eeglab2021.1'); %EEGLAB toolbox
addpath('/Users/claraziane/OneDrive - Universite de Montreal/S2M/projetRAC/dataAnalysis/Toolbox/GED-master'); %For Gaussian filtering

Participants = {'sub-001'; 'sub-002'; 'sub-003'; 'sub-004'; 'sub-005'; 'sub-006'; 'sub-007'; 'sub-008'; 'sub-009'; 'sub-010'; 'sub-011'; 'sub-012'; 'sub-013'; 'sub-014'; 'sub-015'; 'sub-016'; 'sub-017'; 'sub-018'; 'sub-019'; 'sub-020'};
Conditions   = {'Pref'; 'Slow'; 'Fast'; 'uncuedPref'; 'cuedPrefRest'};

% Electrode used for 'best-electrode' analyses
electrode = 'cz';

% Parameters for eigendecomposition
lowFreq  = 1;   % Distance of neighboring frequencies away from stim frequency
highFreq = 1;   % Distance of neighboring frequencies away from stim frequency
sFWHM    = 0.5;   % FWHM of stim frequency
lowFWHM  = 1; % FWHM of the neighboring frequencies
highFWHM = 1; % FWHM of the neighboring frequencies

% Load results previously computed
load([pathResults 'subAll/resultsEEG.mat']); %Comment if running script for the 1st time

eeglab;
for iParticipant = 1:length(Participants)
    participantStr = strcat('SUB', Participants{iParticipant}(end-2:end));
    disp(participantStr)

    % Load data previously extracted for participant
    load([pathData Participants{iParticipant} '/eeg/equalDATA.mat']);

    % Find condition with least amount of trials to homogenize amount of data per condition
    for iCondition = 1:length(Conditions)-1
        dataSize(iCondition) = size(dataAll.([Conditions{iCondition}]),3);
    end
    dataMin = min(dataSize);

    for iCondition = 1:length(Conditions)

        % Create folder in participant's result folder if does not exist
        directoryResults = fullfile(pathResults, Participants{iParticipant}, '/Stable/', Conditions{iCondition}, '/');
        if ~exist(directoryResults, 'dir')
            mkdir(directoryResults)
        end

        % Load data
        data           = dataAll.([Conditions{iCondition}]);
        dataTempLength = dataLength(~isnan(dataLength(:,iCondition)),iCondition);
        if strcmpi(Conditions{iCondition}, 'cuedPrefRest') ~= 1
            stepOnset = stepData.([Conditions{iCondition}]);
        end
        if strcmpi(Conditions{iCondition}, 'uncuedPref') ~= 1
            cueOnset = cueData.([Conditions{iCondition}]);
        end
        freqTemp = freqAll(:,iCondition);
        if size(data,3) > dataMin && strcmpi(Conditions{iCondition}, 'cuedPrefRest') ~= 1
                n2remove = size(data,3) - dataMin;
                trials2remove = round(linspace(1,size(data,3),n2remove));
                data(:,:,trials2remove) = [];
                dataTempLength(trials2remove) = [];
                freqTemp(trials2remove) = [];
                stepOnset(:,trials2remove) = [];
                if strcmpi(Conditions{iCondition}, 'uncuedPref') ~= 1
                    cueOnset(:,trials2remove) = [];
                end
        else
            trials2remove = [];
        end
        nChan = length(chanLocs);

        % FFT Parameters
        fftRes   = ceil(rateEEG/.02); % FFT resolution of .02 Hz
        Hz       = linspace(0,rateEEG,fftRes);

        %% Computing covariiance matrices
        compTime = nan(size(data,2), size(data,3));
        compPhaseHz = nan(size(data,2), size(data,3));
        compFiltered = nan(size(data,2), size(data,3));
        for iBlock = 1:size(data,3)

            if ~isnan(freqTemp(iBlock))
                dataTemp = squeeze(data(:,~isnan(squeeze(data(1,:,iBlock))),iBlock));
                dataFFT(:,:,iBlock) = abs(fft(dataTemp,fftRes,2)/(dataTempLength(iBlock)-1)).^2;

                % S covariance
                sFreq = round(freqTemp(iBlock),2);
                sTemp = filterFGx(dataTemp,rateEEG,sFreq,sFWHM);
                sTemp = bsxfun(@minus,sTemp,mean(sTemp,2));
                sCovariance = (sTemp*sTemp')/(dataTempLength(iBlock)-1);

                % R covariance (low)
                lowTemp = filterFGx(dataTemp,rateEEG,sFreq-lowFreq,lowFWHM);
                lowTemp = bsxfun(@minus,lowTemp,mean(lowTemp,2));
                lowCovariance = (lowTemp*lowTemp')/(dataTempLength(iBlock)-1);

                % R covariance (high)
                highTemp = filterFGx(dataTemp,rateEEG,sFreq+highFreq,highFWHM);
                highTemp = bsxfun(@minus,highTemp,mean(highTemp,2));
                highCovariance = (highTemp*highTemp')/(dataTempLength(iBlock)-1);

                % R covariance (mean of high and low covariance matrices)
                rCovariance    = (highCovariance+lowCovariance)/2;

                % Apply regularization to R (optional)
                regulFactor = .01;
                rCovariance = (1-regulFactor)*rCovariance + regulFactor*mean(eig(rCovariance))*eye(size(rCovariance));

                % Plot covariance martices
%                 clim = [-1 1]*10;
%                 figure(3), clf
%                 subplot(221); imagesc(sCovariance);
%                     axis square; set(gca,'clim',clim); xlabel('Channels'), ylabel('Channels'); colorbar
%                     title('Covariance Matrix S');
%                 subplot(222); imagesc(rCovariance);
%                     axis square; set(gca,'clim',clim); xlabel('Channels'), ylabel('Channels'); colorbar
%                     title('Covariance Matrix R');
%                 subplot(223); imagesc(lowCovariance);
%                     axis square; set(gca,'clim',clim); xlabel('Channels'), ylabel('Channels'); colorbar
%                     title('Covariance Matrix Low (R)');
%                 subplot(224); imagesc(highCovariance);
%                     axis square; set(gca,'clim',clim); xlabel('Channels'), ylabel('Channels'); colorbar
%                     title('Covariance Matrix High (R)');
% %                 saveas(figure(3), [directoryResults '/fig_covarianceMatrices.png']);

                % Generalized eigendecomposition
                [W,Lambdas] = eig(sCovariance, rCovariance);
                [lambdaSorted, lambdaIndex] = sort(diag(Lambdas), 'descend');
                W = bsxfun(@rdivide, W, sqrt(sum(W.^2,1))); % Normalize vectors
                compMaps = sCovariance * W / (W' * sCovariance * W); % Extract components

                % Plot first 5 components
                i = 1;
                figure(4), clf
                subplot(211); plot(lambdaSorted,'ks-','markersize',10,'markerfacecolor','w');
                    xlabel('Component', 'FontSize', 14); ylabel('\lambda', 'FontSize', 14);
                    title('Eigen Values', 'FontSize', 14);
                for iComp = [1:4 length(Lambdas)]

                    % Force Cz to be positive to facilitate across-subject interpretation
                    elecSign = sign(compMaps(elecLoc,lambdaIndex(iComp)));% *(-1) to reverse sign
                    compMaps(:,lambdaIndex(iComp)) = compMaps(:,lambdaIndex(iComp))* elecSign;

                    subplot(2,5,5+i); comp2plot = compMaps(:,lambdaIndex(iComp)); topoplot(comp2plot./max(comp2plot), chanLocs, 'maplimits', [-1 1], 'numcontour',0,'electrodes','on','shading','interp');
                        title(['Component ' num2str(iComp)], 'FontSize', 14)
                    i = i+1;

                end
                colormap jet
% %                 saveas(figure(4), [directoryResults '/fig_ressComponents.png']);

                %% Reconstruct component time series

                comp2Keep(:,iBlock) = input('Which component should be kept ?');%resultsEEG.(participantStr).([Conditions{iCondition}]).compMax(iBlock);
                compMax             = lambdaIndex(comp2Keep(:,iBlock));
                Map(:,iBlock)       = compMaps(:,compMax);

                compTime(1:dataTempLength(iBlock),iBlock) = W(:,compMax)'*squeeze(data(:,~isnan(data(1,:,iBlock)),iBlock));
                compFFT(:,iBlock) = abs(fft(compTime(~isnan(compTime(:,iBlock)),iBlock),fftRes,1)/(dataTempLength(iBlock)-1) ).^2;

                [M, fIndex]     = max(compFFT(:,iBlock));
                freqMax(iBlock) = round(Hz(fIndex),2);
                % If power peak at a frequency due to artifact (i.e., subject 17, condition 5)
%                 if round(freqMax(iBlock)/sFreq) ~= 1
%                     [LIA minLoc] = ismember(sFreq/2, round(Hz,2));
%                     [LIA maxLoc] = ismember(sFreq*2, round(Hz,2));
%                     [M, fIndex]     = max(compFFT(minLoc:maxLoc,iBlock));
%                     freqMax(iBlock) = round(Hz(fIndex+minLoc-1),2);
%                 end
                timeVector      = linspace(1, round(length(compTime)/rateEEG), length(compTime));

                figure(5)
                subplot(221); comp2plot = compMaps(:,compMax); topoplot(comp2plot./max(comp2plot), chanLocs, 'maplimits', [-1 1], 'numcontour',0,'conv','off','electrodes','on','shading','interp'); colorbar;
                    title('Component Topography', 'FontSize', 14);
                subplot(222); plot(timeVector, compTime(:,iBlock));
                    set(gca, 'xlim', [timeVector(1) timeVector(end)]);
                    xlabel({'Time (s)'}, 'FontSize', 14),
                    title('Component Time Series', 'FontSize', 14);
                subplot(2,2,[3:4]); plot(Hz,compFFT(:,iBlock));
                    xlim = [0 25]; set(gca,'xlim',xlim);...
                    xlabel('Frequency (Hz)', 'FontSize', 14) ; ylabel('Power', 'FontSize', 14);
                    legend(['Peak frequency = ' num2str(freqMax(iBlock))], 'FontSize', 14);
                title('Component FFT', 'FontSize', 14);
% %                 saveas(figure(5), [directoryResults '/fig_ressTopo.png']);

                %% Compute SNR spectrum
                elecFFT = dataFFT(elecLoc,:,iBlock);

                [compSNRTemp,elecSNRTemp] = deal(zeros(size(Hz)));
                bins2skip =  5;
                binsNb    = 20+bins2skip;

                % Loop over frequencies to compute SNR
                for iHz = binsNb+1:length(Hz)-binsNb-1
                    numer = compFFT(iHz,iBlock);
                    denom = mean(compFFT([iHz-binsNb:iHz-bins2skip iHz+bins2skip:iHz+binsNb],iBlock));
                    compSNRTemp(iHz) = numer./denom;

                    numer = elecFFT(iHz);
                    denom = mean( elecFFT([iHz-binsNb:iHz-bins2skip iHz+bins2skip:iHz+binsNb]) );
                    elecSNRTemp(iHz) = numer./denom;
                end
                compSNR(:,iBlock) = compSNRTemp;
                elecSNR(:,iBlock) = elecSNRTemp;

%                 figure(6), clf
%                 xlim = [0.5 15];
%                 subplot(2,2,1); comp2plot = compMaps(:,compMax); topoplot(comp2plot./max(comp2plot),chanLocs,'maplimits',[-1 1],'numcontour',0,'electrodes','on','shading','interp');
%                     title([ 'Component at ' num2str(freqMax(iBlock)) ' Hz' ], 'FontSize', 14);
%                 subplot(2,2,2); comp2plot = dataFFT(:,dsearchn(Hz', freqMax(iBlock))); topoplot(comp2plot./max(comp2plot),chanLocs,'maplimits',[-1 1],'numcontour',0,'electrodes','on','emarker2',{find(strcmpi({chanLocs.labels},electrode)) 'o' 'w' 4},'shading','interp');
%                     title([ 'Electrode Power at ' num2str(freqMax(iBlock)) ' Hz' ], 'FontSize', 14);
%                 subplot(2,2,[3:4]); plot(Hz,compSNRTemp,'ro-','linew',1,'markersize',5,'markerface','w'); hold on;
%                     plot(Hz,elecSNRTemp,'ko-','linew',1,'markersize',5,'markerface','w');
%                     set(gca,'xlim',xlim); xlabel('Frequency (Hz)', 'FontSize', 14), ylabel('SNR', 'FontSize', 14); legend({'Component'; electrode}, 'FontSize', 14); clear xlim
% %                 saveas(figure(6), [directoryResults '/fig_ressVSelectrode.png']);

                % Filter data to compute instantaneous frequencies
                compFiltered(1:length(compTime(~isnan(compTime(:,iBlock)),iBlock)),iBlock)  = filterFGx(compTime(~isnan(compTime(:,iBlock)),iBlock)', rateEEG, freqMax(iBlock), sFWHM);

                clear compHilbert compMaps compPhaseAngle...
                    compSNRTemp  elecSNRTemp elecFFT Index highCovariance...
                    lambdaIndex Lambdas lambdaSorted lowCovariance...
                    rCovariance sCovariance W
                close all

            end

        end % End of Block

        %% Compute stability inde

        % Compute Hilbert Transform
        compFiltered = reshape(compFiltered, 1, []);
        compFiltered(isnan(compFiltered)) = [];
        compHilbert = hilbert(compFiltered);

        % Extract phase angles at each step and beat
        compPhase = angle(compHilbert);
%         figure; plot(compPhase); hold on;

        stepPhase = [];
        beatPhase = [];
        timeStartBlock = 1; 
        for iBlock = 1:size(data,3)

            timeEndBlock = timeStartBlock + (dataTempLength(iBlock)-1);
            if strcmpi(Conditions{iCondition}, 'cuedPrefRest') ~= 1
                for iStep = 2:size(stepOnset,1)-1
                    if iStep == 2
                        stepTime  = timeStartBlock + diff(stepOnset(iStep-1:iStep,iBlock));
                    else
                        stepTime  = stepTime + diff(stepOnset(iStep-1:iStep,iBlock)); 
                    end
                    stepPhase = [stepPhase; compPhase(stepTime)];

                    if strcmpi(Conditions{iCondition}, 'uncuedPref') ~= 1
                        [minTime indexTime] = min(abs(cueOnset(:,iBlock)-stepOnset(iStep,iBlock)));
                        minDir = cueOnset(indexTime,iBlock) - stepOnset(iStep,iBlock);
                        if minDir >= 0
                            beatTime  = stepTime+minTime;
                        else
                            beatTime  = stepTime-minTime;
                        end
                        if beatTime > 0
                            beatPhase = [beatPhase; compPhase(beatTime)]; plot([beatTime beatTime], [-3 3], 'r-'); hold on;
                        end
                    end
    
                end
            else
                for iBeat = 1:sum(~isnan(cueOnset(:,end)))
                    if iBeat == 1
                        beatTime = timeStartBlock;
                    else
                        beatTime = timeStartBlock + diff(cueOnset(iBeat-1:iBeat,iBlock));
                    end
                    beatPhase = [beatPhase; compPhase(beatTime)];
                end

            end
            timeStartBlock = timeEndBlock+1;

        end

        % Convert phase angles to Hz
        compPhaseAngle = unwrap(compPhase);
        compPhaseHz = (rateEEG*diff(compPhaseAngle)) / (2*pi);

        % Apply a sliding moving median with a window width of 400 samples
        nOrder = 10;
        orders = linspace(10,400,nOrder)/2;
        orders = round(orders/(1000/rateEEG));

        for iRound = 1
            if iRound == 1
                data2filt = compPhaseHz;
            else
                data2filt = mean(phaseMed);
            end

            phaseMed = zeros(length(orders), length(data2filt));
            for iOrder = 1:nOrder
                for iTime = 1:length(data2filt)
                    phaseTemp = sort(data2filt(max(iTime-orders(iOrder),1):min(iTime+orders(iOrder),length(data2filt)-1)));
                    phaseMed(iOrder,iTime) = phaseTemp(floor(numel(phaseTemp)/2)+1);
                end
            end

        end
        phaseMedFilt = mean(phaseMed);

        % Compute stability index
        timeStartTemp = 1;            
        for iBlock = 1:size(data,3)

            if strcmpi(Conditions{iCondition}, 'cuedPrefRest')
                timeEndTemp = timeStartTemp + (dataTempLength(iBlock)-1);
                stabilityIndex(iBlock) = std(phaseMedFilt(timeStartTemp:timeEndTemp-1));
                timeStartTemp = timeEndTemp+1;                
            else
                timeStart   = timeStartTemp + diff(stepOnset(1:2,iBlock));
                timeEndTemp = timeStartTemp + (dataTempLength(iBlock)-1);
                timeEnd     = timeEndTemp - diff(stepOnset(end-1:end,iBlock));
                stabilityIndex(iBlock) = std(phaseMedFilt(timeStart:timeEnd));
                timeStartTemp = timeEndTemp+1;
            end

        end

        % Convert stability indexes to z-scores to find outliers
        siZ = (stabilityIndex - mean(stabilityIndex)) / std(stabilityIndex);
        siOutliers = abs(siZ) >= 3;
        
        % Remove outliers
        stabilityIndex(siOutliers) = [];

        % Calculate the mean;
        stabilityIndex = mean(stabilityIndex);

        % Plot instantaneous frequencies
        figure(7);
            time = linspace(0, round(length(phaseMedFilt)/rateEEG), length(phaseMedFilt));
%             plot(time, compPhaseHz, 'r--'); hold on;
            plot(time, phaseMedFilt, 'k-'); hold on;
            set(gca, 'xlim', [time(1) time(end)]);
            limY = get(gca, 'ylim');
            plot([1 time(end)], [nanmean(freqTemp) nanmean(freqTemp)], 'color', [0.80,0.80,0.80]); hold on;
            legend({'After moving median smoothing',  'Mean step frequency'}, 'FontSize', 14); %'Before moving median smoothing',
            xlabel({'Time (s)'}, 'FontSize', 14), ylabel({'Frequency (Hz)'}, 'FontSize', 14);
            txt = (['Stability index = ' num2str(mean(stabilityIndex))]); dim = [.2 .5 .3 .3]; annotation('textbox',dim,'String',txt, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FitBoxToText','on', 'FontSize', 14);
            title('Instantaneous frequencies of the extracted component', 'FontSize', 16)
            saveas(figure(7), [directoryResults '/fig_stabilityIndex.png']);

        % Save EEG results
        resultsEEG.(participantStr).([Conditions{iCondition}]).stabilityIndex = stabilityIndex;
        resultsEEG.(participantStr).([Conditions{iCondition}]).compPhase      = compPhase;
        resultsEEG.(participantStr).([Conditions{iCondition}]).stepPhase      = stepPhase;
        resultsEEG.(participantStr).([Conditions{iCondition}]).beatPhase      = beatPhase;
        resultsEEG.(participantStr).([Conditions{iCondition}]).phaseHz        = phaseMedFilt;
        resultsEEG.(participantStr).([Conditions{iCondition}]).compTime       = compTime;
        resultsEEG.(participantStr).([Conditions{iCondition}]).Map            = Map;
        resultsEEG.(participantStr).([Conditions{iCondition}]).chanLocs       = chanLocs;
        resultsEEG.(participantStr).([Conditions{iCondition}]).compSNR        = compSNR;
        resultsEEG.(participantStr).([Conditions{iCondition}]).elecSNR        = elecSNR;
        resultsEEG.(participantStr).([Conditions{iCondition}]).compMax        = comp2Keep;
        resultsEEG.(participantStr).([Conditions{iCondition}]).freqEEG        = rateEEG;
        resultsEEG.(participantStr).([Conditions{iCondition}]).freqMax        = freqMax;
        resultsEEG.(participantStr).([Conditions{iCondition}]).trialsRemoved  = trials2remove;
        resultsEEG.(participantStr).([Conditions{iCondition}]).freqAll        = freqTemp;

        save([pathResults 'subAll/resultsEEG'], 'resultsEEG');

        clear compPhase freqMax compFiltered compFFT compPhaseHz compTime trials2remove dataTempLength dataFFT...
            compSNR elecSNR data phaseMed phaseMedFilt phaseTemp time freqTemp Hz Map comp2Keep timeStartTemp timeEndTemp timeStart timeEnd stepOnset cueOnset stepPhase beatPhase
        close all

    end

    clear dataAll dataLength freqAll elecLoc rateEEG chanLocs

end