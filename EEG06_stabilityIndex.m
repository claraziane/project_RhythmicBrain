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
pathResults = '/Volumes/Seagate/project_rhythmicBrain/Results/';  %Folder where to save results
addpath('/Users/claraziane/Documents/Académique/Informatique/MATLAB/eeglab2021.1'); %EEGLAB toolbox
addpath('/Users/claraziane/OneDrive - Universite de Montreal/S2M/projetRAC/dataAnalysis/Toolbox/GED-master'); %For Gaussian filtering

Participants = {'sub-001'; 'sub-002'; 'sub-003'; 'sub-004'; 'sub-005'; 'sub-006'; 'sub-007'; 'sub-008'; 'sub-009'; 'sub-010'; 'sub-011'; 'sub-012'; 'sub-013'; 'sub-014'; 'sub-015'; 'sub-016'; 'sub-017'; 'sub-018'};
Conditions     = {'Pref'; 'Slow'; 'Fast'; 'uncuedPref'; 'cuedPrefRest'};

% Electrode used for 'best-electrode' analyses
electrode = 'cz';

% Parameters for eigendecomposition
lowFreq  = 1;   % Distance of neighboring frequencies away from stim frequency
highFreq = 1;   % Distance of neighboring frequencies away from stim frequency
sFWHM    = 0.5; % FWHM of stim frequency
lowFWHM  = 0.5; % FWHM of the neighboring frequencies
highFWHM = 0.5; % FWHM of the neighboring frequencies

% Load data and results previously computed
load([pathResults 'subAll/DATA.mat']);
load([pathResults 'subAll/resultsEEG.mat']); %Comment if running script for the 1st time

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for iParticipant = 1:length(Participants)
    participantStr = strcat('SUB', Participants{iParticipant}(end-2:end));

    for iCondition = 1:length(Conditions)

        % Create folder in participant's result folder if does not exist
        directoryResults = fullfile(pathResults, Participants{iParticipant}, '/', Conditions{iCondition}, '/');
        if ~exist(directoryResults, 'dir')
            mkdir(directoryResults)
        end
        
        % Load data
        dataTemp       = dataAll.(participantStr).([Conditions{iCondition}]);
        dataTempLength = dataLength.(participantStr)(~isnan(dataLength.(participantStr)(:,iCondition)),iCondition);
        czLoc          = elecLoc.(participantStr);
        nChan          = length(chanLocs.(participantStr));

        % Keep length of signal corresponding to smallest trial only (avoid transients + homogenize trial length)
        [minLength, indexLength] = min(dataTempLength);
        timeWin = [1 minLength];
        for i = 1:size(dataTemp,3)
            dataTemp(:,end-dataTempLength(i)+1:end,i) = dataTemp(:,~isnan(squeeze(dataTemp(1,:,i))),i);
        end
        dataTemp(:,1:size(dataTemp,2)-minLength,:) = [];

        %% Compute covariance matrices
        % S covariance
        sFreq = round(nanmean(freqAll.(participantStr)(:,iCondition)),2);
        sData = filterFGx(dataTemp,rateEEG.(participantStr),sFreq,sFWHM);
        sData = reshape(sData, nChan,[]);
        sData = bsxfun(@minus,sData,mean(sData,2));
        sCovariance = (sData*sData')/diff(timeWin);

        % R covariance (low)
        lowData = filterFGx(dataTemp,rateEEG.(participantStr),sFreq-lowFreq,lowFWHM);
        lowData = reshape(lowData, nChan,[]);
        lowData = bsxfun(@minus,lowData,mean(lowData,2));
        lowCovariance = (lowData*lowData')/diff(timeWin);

        % R covariance (high)
        highData = filterFGx(dataTemp,rateEEG.(participantStr),sFreq+highFreq,highFWHM);
        highData = reshape(highData, nChan,[]);
        highData = bsxfun(@minus,highData,mean(highData,2));
        highCovariance = (highData*highData')/diff(timeWin);

        % R covariance (mean of high and low covariance matrices)
        rCovariance    = (highCovariance+lowCovariance)/2;

        %% FFT
        fftRes   = ceil(rateEEG.(participantStr)/.02); % FFT resolution of .02 Hz
        Hz       = linspace(0,rateEEG.(participantStr),fftRes);
        dataFFT  = mean(abs(fft(dataTemp,fftRes,2)/diff(timeWin)).^2,3);

        % Apply regularization to R (optional)
        regulFactor = .01;
        rCovariance = (1-regulFactor)*rCovariance + regulFactor*mean(eig(rCovariance))*eye(size(rCovariance));

        % Plot covariance martices
        clim = [-1 1]*1000;
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
            elecSign = sign(compMaps(czLoc,lambdaIndex(iComp)));% *(-1) to reverse sign
            compMaps(:,lambdaIndex(iComp)) = compMaps(:,lambdaIndex(iComp))* elecSign;

            subplot(2,5,5+i); comp2plot = compMaps(:,lambdaIndex(iComp)); topoplot(comp2plot./max(comp2plot), chanLocs.(participantStr), 'maplimits', [-1 1], 'numcontour',0,'electrodes','off','shading','interp');
            title(['Component ' num2str(iComp)], 'FontSize', 14)
            i = i+1;

        end
        colormap jet
        saveas(figure(4), [directoryResults '/fig_ressComponents.png']);

        %% Reconstruct component time series

        comp2Keep   = 1;%input('Which component should be kept ?');
        compMax     = lambdaIndex(comp2Keep);

        compTime = nan(size(dataTemp,2), size(dataTemp,3));
        for i = 1:size(dataTemp,3)
            compTime(:,i) = W(:,compMax)'*squeeze(dataTemp(:,:,i));
        end

        compFFT = mean(abs(fft(compTime,fftRes,1)/diff(timeWin) ).^2,2);
        [M, fIndex] = max(compFFT);
        timeVector = linspace(1, round(length(compTime)/rateEEG.(participantStr)), length(compTime));

        figure(5)
        subplot(221); comp2plot = compMaps(:,compMax); topoplot(comp2plot./max(comp2plot), chanLocs.(participantStr), 'maplimits', [-1 1], 'numcontour',0,'conv','off','electrodes','on','shading','interp'); colorbar;
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
        elecFFT    = dataFFT(czLoc,:,:);

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
        subplot(2,2,1); comp2plot = compMaps(:,compMax); topoplot(comp2plot./max(comp2plot),chanLocs.(participantStr),'maplimits',[-1 1],'numcontour',0,'electrodes','on','shading','interp');
        title([ 'Component for ' num2str(sFreq) ' Hz' ], 'FontSize', 14);
        subplot(2,2,2); comp2plot = dataFFT(:,dsearchn(Hz', sFreq)); topoplot(comp2plot./max(comp2plot),chanLocs.(participantStr),'maplimits',[-1 1],'numcontour',0,'electrodes','on','emarker2',{find(strcmpi({chanLocs.(participantStr).labels},electrode)) 'o' 'w' 4},'shading','interp');
        title([ 'Electrode Power at ' num2str(sFreq) ' Hz' ], 'FontSize', 14);
        subplot(2,2,[3:4]); plot(Hz,compSNR,'ro-','linew',1,'markersize',5,'markerface','w'); hold on;
        plot(Hz,elecSNR,'ko-','linew',1,'markersize',5,'markerface','w');
        set(gca,'xlim',xlim); xlabel('Frequency (Hz)', 'FontSize', 14), ylabel('SNR', 'FontSize', 14); legend({'Component'; electrode}, 'FontSize', 14); clear xlim
        saveas(figure(6), [directoryResults '/fig_ressVSelectrode.png']);

        %% Compute stability index

        % Gaussian filering extracted component
        compTimeConcat = reshape(compTime, 1,[]);
        compFiltered   = filterFGx(compTimeConcat, rateEEG.(participantStr), sFreq, 1);

        % Compute Hilbert Transform
        compHilbert = hilbert(compFiltered);

        % Extract phase angles
        compPhaseAngle = unwrap(angle(compHilbert));
        compPhaseHz    = (rateEEG.(participantStr)*diff(compPhaseAngle)) / (2*pi);

        % Apply a sliding moving median with a window width of 400 samples
        nOrder = 10;
        orders = linspace(10,400,nOrder)/2;
        orders = round(orders/(1000/rateEEG.(participantStr)));
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
        time = linspace(1, round(length(phaseMedFilt)/rateEEG.(participantStr)), length(phaseMedFilt));
        plot(time, compPhaseHz, 'r--'); hold on;
        plot(time, phaseMedFilt, 'k-'); hold on;
        set(gca, 'xlim', [time(1) time(end)]);
        limY = get(gca, 'ylim');
        plot([1 time(end)], [sFreq sFreq], 'color', [0.80,0.80,0.80]); hold on;
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
        resultsEEG.(participantStr).([Conditions{iCondition}]).chanLocs       = chanLocs.(participantStr);
        resultsEEG.(participantStr).([Conditions{iCondition}]).compSNR        = compSNR;
        resultsEEG.(participantStr).([Conditions{iCondition}]).elecSNR        = elecSNR;
        resultsEEG.(participantStr).([Conditions{iCondition}]).compMax        = comp2Keep;
        resultsEEG.(participantStr).([Conditions{iCondition}]).freqEEG        = rateEEG.(participantStr);
        resultsEEG.(participantStr).([Conditions{iCondition}]).freqStep       = sFreq;
        save([pathResults 'subAll/resultsEEG'], 'resultsEEG');

        clear comp2plot compFFT compFiltered compHilbert compMaps compPhaseAngle...
            compPhaseHz compSNR compTime compTimeConcat...
            dataTemp dataTempLength dataFFT elecSNR elecFFT Index...
            highCovariance highData Hz lambdaIndex Lambdas lambdaSorted lowCovariance lowData...
            phaseMed phaseMedFilt phaseTemp rCovariance sCovariance sData time timeVector W
        close all

    end

    clear czLoc

end