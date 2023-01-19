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

Participants = {'sub-001'; 'sub-002'; 'sub-003'; 'sub-004'; 'sub-005'; 'sub-006'; 'sub-007'; 'sub-008'; 'sub-009'; 'sub-010'; 'sub-011'; 'sub-012'; 'sub-013'; 'sub-014'; 'sub-015'; 'sub-016'; 'sub-017'; 'sub-018'};
Conditions   = {'Pref'; 'Slow'; 'Fast'; 'uncuedPref'; 'cuedPrefRest'};

% Electrode used for 'best-electrode' analyses
electrode = 'cz';

% Parameters for eigendecomposition
lowFreq  = 1;   % Distance of neighboring frequencies away from stim frequency
highFreq = 1;   % Distance of neighboring frequencies away from stim frequency
sFWHM    = 0.5; % FWHM of stim frequency
lowFWHM  = 0.5; % FWHM of the neighboring frequencies
highFWHM = 0.5; % FWHM of the neighboring frequencies

% Load results previously computed
% load([pathResults 'subAll/DATA.mat']);
load([pathResults 'subAll/resultsEEG.mat']); %Comment if running script for the 1st time

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for iParticipant = 1:length(Participants)
    participantStr = strcat('SUB', Participants{iParticipant}(end-2:end));

    % Load data previously extracted for participant
    load([pathData Participants{iParticipant} '/eeg/DATA.mat']);

    for iCondition = 1:length(Conditions)

        % Create folder in participant's result folder if does not exist
        directoryResults = fullfile(pathResults, Participants{iParticipant}, '/', Conditions{iCondition}, '/');
        if ~exist(directoryResults, 'dir')
            mkdir(directoryResults)
        end
        
        % Load data
        data           = dataAll.([Conditions{iCondition}]);
        dataTempLength = dataLength(~isnan(dataLength(:,iCondition)),iCondition);
        nChan          = length(chanLocs);

        % FFT Parameters
        fftRes   = ceil(rateEEG/.02); % FFT resolution of .02 Hz
        Hz       = linspace(0,rateEEG,fftRes);
        
        %% Computing covariiance matrices
        sData    = nan(nChan, max(dataTempLength), size(data,3)); %Preallocating matrix
        lowData  = nan(nChan, max(dataTempLength), size(data,3)); %Preallocating matrix
        highData = nan(nChan, max(dataTempLength), size(data,3)); %Preallocating matrix     
        for iBlock = 1:size(data,3)

            if ~isnan(freqAll(iBlock,iCondition))
                dataTemp = squeeze(data(:,~isnan(squeeze(data(1,:,iBlock))),iBlock));
                dataFFT(:,:,iBlock) = abs(fft(dataTemp,fftRes,2)/(dataTempLength(iBlock)-1)).^2;

                % S covariance
                sFreq = round(freqAll(iBlock,iCondition),2);
                sTemp = filterFGx(dataTemp,rateEEG,sFreq,sFWHM);
                sTemp = bsxfun(@minus,sTemp,mean(sTemp,2));
                sCovariance(:,:,iBlock) = (sTemp*sTemp')/(dataTempLength(iBlock)-1);
                sData(:,1:length(sTemp),iBlock) = sTemp;

                % R covariance (low)
                lowTemp = filterFGx(dataTemp,rateEEG,sFreq-lowFreq,lowFWHM);
                lowTemp = bsxfun(@minus,lowTemp,mean(lowTemp,2));
                lowCovariance(:,:,iBlock) = (lowTemp*lowTemp')/(dataTempLength(iBlock)-1);
                lowData(:,1:length(lowTemp),iBlock) = lowTemp;

                % R covariance (high)
                highTemp = filterFGx(dataTemp,rateEEG,sFreq+highFreq,highFWHM);
                highTemp = bsxfun(@minus,highTemp,mean(highTemp,2));
                highCovariance(:,:,iBlock) = (highTemp*highTemp')/(dataTempLength(iBlock)-1);
                highData(:,1:length(highTemp),iBlock) = highTemp;

                clear dataTemp sTemp lowTemp highTemp
            end

        end
        dataFFT = mean(dataFFT,3);



        sCovarianceMean    = mean(sCovariance,3);
        lowCovarianceMean  = mean(lowCovariance,3);
        highCovarianceMean = mean(highCovariance,3);

        % Compute Euclidean distance
        sCovarianceDistance    = zeros(size(sCovariance,3),1);
        lowCovarianceDistance  = zeros(size(lowCovariance,3),1);
        highCovarianceDistance = zeros(size(highCovariance,3),1);
        for iBlock = 1:size(sCovariance,3)
            sCovarianceTemp = squeeze(sCovariance(:,:,iBlock));
            sCovarianceDistance(iBlock) = sqrt(sum((sCovarianceTemp(:) - sCovarianceMean(:)).^2));

            lowCovarianceTemp = squeeze(lowCovariance(:,:,iBlock));
            lowCovarianceDistance(iBlock) = sqrt(sum((lowCovarianceTemp(:) - lowCovarianceMean(:)).^2));

            highCovarianceTemp = squeeze(highCovariance(:,:,iBlock));
            highCovarianceDistance(iBlock) = sqrt(sum((highCovarianceTemp(:) - highCovarianceMean(:)).^2));        
        end
        sCovarianceZ    = (sCovarianceDistance-mean(sCovarianceDistance))       / std(sCovarianceDistance);    %Convert to Z scores
        lowCovarianceZ  = (lowCovarianceDistance-mean(lowCovarianceDistance))   / std(lowCovarianceDistance);  %Convert to Z scores
        highCovarianceZ = (highCovarianceDistance-mean(highCovarianceDistance)) / std(highCovarianceDistance); %Convert to Z scores

%         % Plot covariance distances
%         figure, clf
%         subplot(3,1,1); plot(sCovarianceZ,'ks-','linew',2,'markerfacecolor','w','markersize',12);
%         xlabel('Block'), ylabel('Z_{distance}'); title('Z-scored covariance S distances')
%         subplot(3,1,2); plot(lowCovarianceZ,'ks-','linew',2,'markerfacecolor','w','markersize',12);
%         xlabel('Block'), ylabel('Z_{distance}'); title('Z-scored covariance R (low) distances')
%         subplot(3,1,3); plot(highCovarianceZ,'ks-','linew',2,'markerfacecolor','w','markersize',12);
%         xlabel('Block'), ylabel('Z_{distance}'); title('Z-scored covariance R (high) distances')

        % Pick threshold to reject covariance matrice
        Threshold = 2.31; %corresponds to p < 0.01

        % Identify blocks that exceed the threshold
        sCovarianceReject    = sCovarianceZ    > Threshold;
        lowCovarianceReject  = lowCovarianceZ  > Threshold;
        highCovarianceReject = highCovarianceZ > Threshold;

        % Remove blocks from covariance matrices and recompute grand average
        sCovariance(:,:,sCovarianceReject) = [];
        sCovariance = mean(sCovariance,3);

        lowCovariance(:,:,lowCovarianceReject) = [];
        lowCovariance = mean(lowCovariance,3);

        highCovariance(:,:,highCovarianceReject) = [];
        highCovariance = mean(highCovariance,3);

        % R covariance (mean of high and low covariance matrices)
        rCovariance    = (highCovariance+lowCovariance)/2;

        % Apply regularization to R (optional)
        regulFactor = .01;
        rCovariance = (1-regulFactor)*rCovariance + regulFactor*mean(eig(rCovariance))*eye(size(rCovariance));        

        % Plot covariance martices
        clim = [-1 1]*10;
        figure(3), clf
        subplot(221); imagesc(sCovariance);
        axis square; set(gca,'clim',clim); xlabel('Channels'), ylabel('Channels'); colorbar
        title('Covariance Matrix S');
        subplot(222); imagesc(rCovariance);
        axis square; set(gca,'clim',clim); xlabel('Channels'), ylabel('Channels'); colorbar
        title('Covariance Matrix R');
        subplot(223); imagesc(lowCovariance);
        axis square; set(gca,'clim',clim); xlabel('Channels'), ylabel('Channels'); colorbar
        title('Covariance Matrix Low (R)');
        subplot(224); imagesc(highCovariance);
        axis square; set(gca,'clim',clim); xlabel('Channels'), ylabel('Channels'); colorbar
        title('Covariance Matrix High (R)');
        saveas(figure(3), [directoryResults '/fig_covarianceMatrices.png']);

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

            subplot(2,5,5+i); comp2plot = compMaps(:,lambdaIndex(iComp)); topoplot(comp2plot./max(comp2plot), chanLocs, 'maplimits', [-1 1], 'numcontour',0,'electrodes','off','shading','interp');
            title(['Component ' num2str(iComp)], 'FontSize', 14)
            i = i+1;

        end
        colormap jet
        saveas(figure(4), [directoryResults '/fig_ressComponents.png']);

        %% Reconstruct component time series

        comp2Keep   = 1;%input('Which component should be kept ?');
        compMax     = lambdaIndex(comp2Keep);

        compTime = nan(size(data,2), size(data,3));
        for iBlock = 1:size(data,3)
            compTime(1:dataTempLength(iBlock),iBlock) = W(:,compMax)'*squeeze(data(:,~isnan(data(1,:,iBlock)),iBlock));
            compFFT(:,iBlock)  = abs(fft(compTime(~isnan(compTime(:,iBlock)),iBlock),fftRes,1)/(dataTempLength(iBlock)-1) ).^2;
        end

        compFFT = mean(compFFT,2);
        [M, fIndex] = max(compFFT);
        timeVector = linspace(1, round(length(compTime)/rateEEG), length(compTime));

        figure(5)
        subplot(221); comp2plot = compMaps(:,compMax); topoplot(comp2plot./max(comp2plot), chanLocs, 'maplimits', [-1 1], 'numcontour',0,'conv','off','electrodes','on','shading','interp'); colorbar;
        title('Component Topography', 'FontSize', 14);
        subplot(222); plot(timeVector, compTime);
        set(gca, 'xlim', [timeVector(1) timeVector(end)]);
        xlabel({'Time (s)'}, 'FontSize', 14),
        title('Component Time Series', 'FontSize', 14);
        subplot(2,2,[3:4]); plot(Hz,compFFT);
        xlim = [0 25]; set(gca,'xlim',xlim);...
            xlabel('Frequency (Hz)', 'FontSize', 14) ; ylabel('Power', 'FontSize', 14);
        legend(['Peak frequency = ' num2str(round(Hz(fIndex),2))], 'FontSize', 14);
        title('Component FFT', 'FontSize', 14);
        saveas(figure(5), [directoryResults '/fig_ressTopo.png']);

        %% Compute SNR spectrum
        elecFFT    = dataFFT(elecLoc,:,:);

        [compSNR,elecSNR] = deal(zeros(size(Hz)));
        bins2skip =  5;
        binsNb    = 20+bins2skip;

        % Loop over frequencies to compute SNR
        for iHz = binsNb+1:length(Hz)-binsNb-1
            numer = compFFT(iHz);
            denom = mean(compFFT([iHz-binsNb:iHz-bins2skip iHz+bins2skip:iHz+binsNb]));
            compSNR(iHz) = numer./denom;

            numer = elecFFT(iHz);
            denom = mean( elecFFT([iHz-binsNb:iHz-bins2skip iHz+bins2skip:iHz+binsNb]) );
            elecSNR(iHz) = numer./denom;
        end

        figure(6), clf
        xlim = [0.5 15];
        subplot(2,2,1); comp2plot = compMaps(:,compMax); topoplot(comp2plot./max(comp2plot),chanLocs,'maplimits',[-1 1],'numcontour',0,'electrodes','on','shading','interp');
        title([ 'Component at ' num2str(round(Hz(fIndex),2)) ' Hz' ], 'FontSize', 14);
        subplot(2,2,2); comp2plot = dataFFT(:,dsearchn(Hz', round(Hz(fIndex),2))); topoplot(comp2plot./max(comp2plot),chanLocs,'maplimits',[-1 1],'numcontour',0,'electrodes','on','emarker2',{find(strcmpi({chanLocs.labels},electrode)) 'o' 'w' 4},'shading','interp');
        title([ 'Electrode Power at ' num2str(round(Hz(fIndex),2)) ' Hz' ], 'FontSize', 14);
        subplot(2,2,[3:4]); plot(Hz,compSNR,'ro-','linew',1,'markersize',5,'markerface','w'); hold on;
        plot(Hz,elecSNR,'ko-','linew',1,'markersize',5,'markerface','w');
        set(gca,'xlim',xlim); xlabel('Frequency (Hz)', 'FontSize', 14), ylabel('SNR', 'FontSize', 14); legend({'Component'; electrode}, 'FontSize', 14); clear xlim
        saveas(figure(6), [directoryResults '/fig_ressVSelectrode.png']);

        %% Compute stability index

        % Gaussian filering extracted component
        compTimeConcat = reshape(compTime, 1,[]);
        compTimeConcat(isnan(compTimeConcat)) = [];
        compFiltered   = filterFGx(compTimeConcat, rateEEG, round(Hz(fIndex),2), 1);

        % Compute Hilbert Transform
        compHilbert = hilbert(compFiltered);

        % Extract phase angles
        compPhaseAngle = unwrap(angle(compHilbert));
        compPhaseHz    = (rateEEG*diff(compPhaseAngle)) / (2*pi);

        % Apply a sliding moving median with a window width of 400 samples
        nOrder = 10;
        orders = linspace(10,400,nOrder)/2;
        orders = round(orders/(1000/rateEEG));
        phaseMed = zeros(length(orders), length(compPhaseHz));
        for iOrder = 1:nOrder
            for iTime = 1:length(compPhaseHz)
                phaseTemp = sort(compPhaseHz(max(iTime-orders(iOrder),1):min(iTime+orders(iOrder),length(compPhaseHz)-1)));
                phaseMed(iOrder,iTime) = phaseTemp(floor(numel(phaseTemp)/2)+1);
            end
        end
        phaseMedFilt = mean(phaseMed);

        % Compute stability index
        stabilityIndex = std(phaseMedFilt);

        % Plot instantaneous frequencies
        figure(7);
        time = linspace(1, round(length(phaseMedFilt)/rateEEG), length(phaseMedFilt));
        plot(time, compPhaseHz, 'r--'); hold on;
        plot(time, phaseMedFilt, 'k-'); hold on;
        set(gca, 'xlim', [time(1) time(end)]);
        limY = get(gca, 'ylim');
        plot([1 time(end)], [round(Hz(fIndex),2) round(Hz(fIndex),2)], 'color', [0.80,0.80,0.80]); hold on;
        legend({'Before moving median smoothing', 'After moving median smoothing',  'Mean step frequency'}, 'FontSize', 14);
        xlabel({'Time (s)'}, 'FontSize', 14), ylabel({'Frequency (Hz)'}, 'FontSize', 14);
        txt = (['Stability index = ' num2str(mean(stabilityIndex))]); dim = [.2 .5 .3 .3]; annotation('textbox',dim,'String',txt, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FitBoxToText','on', 'FontSize', 14);
        title('Instantaneous frequencies of the extracted component', 'FontSize', 16)
        saveas(figure(7), [directoryResults '/fig_stabilityIndex.png']);

        % Save EEG results
        resultsEEG.(participantStr).([Conditions{iCondition}]).stabilityIndex = stabilityIndex;
        resultsEEG.(participantStr).([Conditions{iCondition}]).phaseAngle     = compPhaseAngle;
        resultsEEG.(participantStr).([Conditions{iCondition}]).phaseHz        = phaseMedFilt;
        resultsEEG.(participantStr).([Conditions{iCondition}]).compTime       = compTime;
        resultsEEG.(participantStr).([Conditions{iCondition}]).Map            = compMaps(:,compMax);
        resultsEEG.(participantStr).([Conditions{iCondition}]).chanLocs       = chanLocs;
        resultsEEG.(participantStr).([Conditions{iCondition}]).compSNR        = compSNR;
        resultsEEG.(participantStr).([Conditions{iCondition}]).elecSNR        = elecSNR;
        resultsEEG.(participantStr).([Conditions{iCondition}]).compMax        = comp2Keep;
        resultsEEG.(participantStr).([Conditions{iCondition}]).freqEEG        = rateEEG;
        resultsEEG.(participantStr).([Conditions{iCondition}]).freqStep       = round(Hz(fIndex),2);
        save([pathResults 'subAll/resultsEEG'], 'resultsEEG');

        clear comp2plot compFFT compFiltered compHilbert compMaps compPhaseAngle...
            highCovarianceDistance highCovarianceMean highCovarianceReject highCovarianceTemp highCovarianceZ...
            lowCovarianceDistance lowCovarianceMean lowCovarianceReject lowCovarianceTemp lowCovarianceZ...
            sCovarianceDistance sCovarianceMean sCovarianceReject sCovarianceTemp sCovarianceZ...
            compPhaseHz compSNR compTime compTimeConcat...
            data dataTempLength dataFFT elecSNR elecFFT Index...
            highCovariance highData Hz lambdaIndex Lambdas lambdaSorted lowCovariance lowData...
            phaseMed phaseMedFilt phaseTemp rCovariance sCovariance sData time timeVector W
        close all

    end

    clear dataAll dataLength freqAll elecLoc rateEEG chanLocs

end