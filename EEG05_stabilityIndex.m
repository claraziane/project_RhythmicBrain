%% Calculating the stability index
% 1. Compute RESS component :
%   a. sequence EEG in trials according to beat onsets (-100:500 ms)
%   b. compute covariance matrices S (from narrow-filtered data) and R (from broad-band signal)
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
pathResults = '/Volumes/Seagate/project_rhythmicBrain/Results/';
addpath('/Users/claraziane/Documents/Académique/Informatique/MATLAB/eeglab2021.1');
addpath('/Users/claraziane/OneDrive - Universite de Montreal/S2M/projetRAC/dataAnalysis/Toolbox/GED-master');
addpath '/Volumes/10.89.24.15/Projet_RAC/dataAnalysis/Scripts/Clara/EEG/Functions'

Participants = {'sub-001'; 'sub-002'; 'sub-003'; 'sub-004'; 'sub-005'; 'sub-006'; 'sub-007'; 'sub-008'; 'sub-009'; 'sub-010'; 'sub-011'; 'sub-012'; 'sub-013'; 'sub-014'; 'sub-015'; 'sub-016'; 'sub-017'; 'sub-018'};
Blocks       = {'run-01'; 'run-02'; 'run-03'; 'run-04'; 'run-05'; 'run-06'; 'run-07'; 'run-08'; 'run-09'; 'run-10'; 'run-11'; 'run-12'; 'run-13'; 'run-14'; 'run-15'; 'run-16'; 'run-17'; 'run-18'; 'run-19'};
speedStr     = {'Pref'; 'Slow'; 'Fast'; 'uncuedPref'; 'cuedPrefRest'};
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

% Electrode used for 'best-electrode' analyses
electrode = 'cz';

% Parameters for eigendecomposition
lowFreq  = 1;   % Distance of neighboring frequencies away from stim frequency
highFreq = 1;   % Distance of neighboring frequencies away from stim frequency
sFWHM    = 0.5; % FWHM of stim frequency
lowFWHM  = 0.5; % FWHM of the neighboring frequencies
highFWHM = 0.5; % FWHM of the neighboring frequencies

load([pathResults 'subAll/resultsEEG.mat']);

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for iParticipant = 2:4%length(Participants)
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

            % Keep only EEG channels
            chanEEG = strcmpi('EEG',{EEG.chanlocs.type});
            nChan = sum(chanEEG);

            % Check that the data  matrix is full rank
            if rank(data) ~= nChan
                warning('Data matrix is not full rank')
            end

            % Find Cz in channels
            elecLoc = strcmpi(electrode,{EEG.chanlocs.labels});

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

                for iSpeed = 1:length(speedStr)-1 % Pref: iSpeed = 1; Slow: iSpeed = 2; Fast: iSpeed = 3; Uncued Pref: iSpeed = 4;

                    if iBlock == 1
                        i = 1; %Total number of trials (all blocks)
                    else
                        i = size((dataAll.([speedStr{iSpeed}])),3) + 1;
                    end

                    for iTrial = 1:length(eventLatency(~isnan(eventLatency(iSpeed,:))))

                        if ~isnan(eventLatency(iSpeed,iTrial))

                            % Find time periods of cued gait at preferred cadence
                            if iSpeed == 1
                                if ~isnan(eventCues(iSpeed, iTrial))
                                    timeWin(iSpeed,:) = [eventLatency(iSpeed,iTrial) eventCues(iSpeed, iTrial)-1];
                                elseif isnan(eventCues(iSpeed, iTrial)) && ismember(EEG.urevent(end).type, [101 151 201 251])
                                    timeWin(iSpeed,:) = [eventLatency(iSpeed,iTrial) EEG.urevent(end).latency];
                                end

                            % Find time periods of cued gait at slow and fast cadences
                            else
                                target = eventCues(iSpeed,:) - eventLatency(iSpeed,iTrial+1);
                                target = target(target<0);
                                if ~isempty(target)
                                    timeWin(iSpeed,:) = [eventLatency(iSpeed,iTrial) eventCues(iSpeed, length(target))];
                                else
                                    Index = ~isnan(eventCues(iSpeed,:));
                                    timeWin(iSpeed,:) = [eventLatency(iSpeed,iTrial)  eventCues(iSpeed,sum(Index))];
                                end
                            end

                            if i == 1
                                dataAll.([speedStr{iSpeed}])(:,:,i) = nan(nChan,diff(timeWin(iSpeed,:))+1);
                            else
                                timeDiff = diff(timeWin(iSpeed,:))+1 - size((dataAll.([speedStr{iSpeed}])),2);

                                if timeDiff > 0
                                    dataAll.([speedStr{iSpeed}])(:,size((dataAll.([speedStr{iSpeed}])),2)+1:size((dataAll.([speedStr{iSpeed}])),2)+timeDiff,1:i-1) = nan;
                                    dataAll.([speedStr{iSpeed}])(:,:,i) = nan(nChan,diff(timeWin(iSpeed,:))+1);
                                elseif timeDiff <= 0
                                    dataAll.([speedStr{iSpeed}])(:,:,i) = nan(nChan,size((dataAll.([speedStr{iSpeed}])),2));
                                end

                            end
                            dataAll.([speedStr{iSpeed}])(:,1:diff(timeWin(iSpeed,:))+1,i) = data(:,timeWin(iSpeed,1):timeWin(iSpeed,2));
                            dataLength(i,iSpeed) = diff(timeWin(iSpeed,:))+1;

                            % Compute step freq
                            stepPrefCued   = [];
                            stepSlowCued   = [];
                            stepFastCued   = [];
                            stepPrefUncued = [];
                            for iLatency = 1:length(EEG.urevent)
                                if EEG.urevent(iLatency).latency >= timeWin(iSpeed,1) && EEG.urevent(iLatency).latency <= timeWin(iSpeed,2)
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

                            if iSpeed == 1
                                stepPrefCued = sort(stepPrefCued);
                                Freq(i,iSpeed) = mean(diff(stepPrefCued));
                            elseif iSpeed == 2
                                stepSlowCued = sort(stepSlowCued);
                                Freq(i,iSpeed) = mean(diff(stepSlowCued));
                            elseif iSpeed == 3
                                stepFastCued = sort(stepFastCued);
                                Freq(i,iSpeed) = mean(diff(stepFastCued));
                            elseif iSpeed == 4
                                stepPrefUncued = sort(stepPrefUncued);
                                Freq(i,iSpeed) = mean(diff(stepPrefUncued));
                            end
                            Freq(i,iSpeed) = round(EEG.srate/Freq(i,iSpeed),1);

                            i = i + 1;
                            clear stepPrefCued stepSlowCued stepFastCued target

                        end

                    end

                end

            % If there is no steps in trial (ie., standing condition)
            elseif sum(ismember([EEG.urevent(1:end).type], Event(Event == 171))) ~= 0
                iSpeed = 5;

                restPrefCued = [];
                if sum(ismember([EEG.urevent(1:end).type], [171 271])) == length(EEG.urevent)
                    for iLatency = 1:length(EEG.urevent)
                        restPrefCued = [restPrefCued; EEG.urevent(iLatency).latency];
                    end
                    
                    % Finding time periods of standing data
                    timeWin(iSpeed,1) = restPrefCued(1);
                    timeWin(iSpeed,2) = restPrefCued(end);

                    % Extracting rest data
                    if k == 1
                        dataAll.([speedStr{iSpeed}])(:,:,k) = nan(nChan,diff(timeWin(iSpeed,:))+1);
                    else
                        
                        timeDiff = diff(timeWin(iSpeed,:))+1 - size((dataAll.([speedStr{iSpeed}])),2);
                        if timeDiff > 0
                            dataAll.([speedStr{iSpeed}])(:,size((dataAll.([speedStr{iSpeed}])),2)+1:size((dataAll.([speedStr{iSpeed}])),2)+timeDiff,1:k-1) = nan;
                            dataAll.([speedStr{iSpeed}])(:,:,k) = nan(nChan,diff(timeWin(iSpeed,:))+1);
                        elseif timeDiff <= 0
                            dataAll.([speedStr{iSpeed}])(:,:,k) = nan(nChan,size((dataAll.([speedStr{iSpeed}])),2));
                        end

                    end
                    dataAll.([speedStr{iSpeed}])(:,1:diff(timeWin(iSpeed,:))+1,k) = data(:,timeWin(iSpeed,1):timeWin(iSpeed,2));
                    dataLength(k,iSpeed) = diff(timeWin(iSpeed,:))+1;

                    Freq(k,iSpeed) = mean(diff(restPrefCued));
                    Freq(k,iSpeed) = round(EEG.srate/Freq(k,iSpeed),1);
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

    dataLength(dataLength == 0) = nan;
    Freq(Freq == 0) = nan;

    for iSpeed = 1:length(speedStr) % Pref: iSpeed = 1; Slow: iSpeed = 2; Fast: iSpeed = 3; Uncued Pref: iSpeed = 4; Cued Rest: iSpeed = 5

        % Create folder in participant's result folder if does not exist
        directoryResults = fullfile(pathResults, Participants{iParticipant}, '/', speedStr{iSpeed}, '/');
        if ~exist(directoryResults, 'dir')
            mkdir(directoryResults)
        end

        dataTemp       = dataAll.([speedStr{iSpeed}]);
        dataTempLength = dataLength(~isnan(dataLength(:,iSpeed)),iSpeed);

        % Keep length of signal corresponding to smallest trial only (avoid transients + homogenize trial length)
        [minLength, indexLength] = min(dataTempLength);
        timeWin = [1 minLength];
        for i = 1:size(dataTemp,3)
            dataTemp(:,end-dataTempLength(i)+1:end,i) = dataTemp(:,~isnan(squeeze(dataTemp(1,:,i))),i);
        end
        dataTemp(:,1:size(dataTemp,2)-minLength,:) = [];

        %% Compute covariance matrices
        % S covariance
        sFreq = round(nanmean(Freq(:,iSpeed)),1);
        sData = filterFGx(dataTemp,EEG.srate,sFreq,sFWHM);
        sData = reshape(sData, nChan,[]);
        sData = bsxfun(@minus,sData,mean(sData,2));
        sCovariance = (sData*sData')/diff(timeWin);

        % R covariance (low)
        lowData = filterFGx(dataTemp,EEG.srate,sFreq-lowFreq,lowFWHM);
        lowData = reshape(lowData, nChan,[]);
        lowData = bsxfun(@minus,lowData,mean(lowData,2));
        lowCovariance = (lowData*lowData')/diff(timeWin);

        % R covariance (high)
        highData = filterFGx(dataTemp,EEG.srate,sFreq+highFreq,highFWHM);
        highData = reshape(highData, nChan,[]);
        highData = bsxfun(@minus,highData,mean(highData,2));
        highCovariance = (highData*highData')/diff(timeWin);

        % R covariance (mean of high and low covariance matrices)
        rCovariance    = (highCovariance+lowCovariance)/2;

        %% FFT
        fftRes   = ceil(EEG.srate/.02); % FFT resolution of .02 Hz
        Hz       = linspace(0,EEG.srate,fftRes);
        dataFFT  = mean(abs(fft(dataTemp,fftRes,2)/diff(timeWin)).^2,3);

        %     % Apply regularization to R (optional)
        %     regulFactor = .01;
        %     rCovariance = (1-regulFactor)*rCovariance + regulFactor*mean(eig(rCovariance))*eye(size(rCovariance));

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
            elecSign = sign(compMaps(elecLoc,lambdaIndex(iComp)));% *(-1) to reverse sign
            compMaps(:,lambdaIndex(iComp)) = compMaps(:,lambdaIndex(iComp))* elecSign;

            subplot(2,5,5+i); comp2plot = compMaps(:,lambdaIndex(iComp)); topoplot(comp2plot./max(comp2plot), EEG.chanlocs, 'maplimits', [-1 1], 'numcontour',0,'electrodes','off','shading','interp');
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
        timeVector = linspace(1, round(length(compTime)/EEG.srate), length(compTime));

        figure(5)
        subplot(221); comp2plot = compMaps(:,compMax); topoplot(comp2plot./max(comp2plot), EEG.chanlocs, 'maplimits', [-1 1], 'numcontour',0,'conv','off','electrodes','on','shading','interp'); colorbar;
        title('Component Topography', 'FontSize', 14);
        subplot(222); plot(timeVector, compTime);
        set(gca, 'xlim', [timeVector(1) timeVector(end)]);
        xlabel({'Time (s)'}, 'FontSize', 14),
        title('Component Time Series', 'FontSize', 14);
        subplot(2,2,[3:4]); plot(Hz,compFFT);
        xlim = [0 25]; set(gca,'xlim',xlim);...
            xlabel('Frequency (Hz)', 'FontSize', 14) ; ylabel('Power', 'FontSize', 14);
        legend(['Peak frequency = ' num2str(round(Hz(fIndex),1))], 'FontSize', 14);
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
        subplot(2,2,1); comp2plot = compMaps(:,compMax); topoplot(comp2plot./max(comp2plot),EEG.chanlocs,'maplimits',[-1 1],'numcontour',0,'electrodes','on','shading','interp');
        title([ 'Component for ' num2str(sFreq) ' Hz' ], 'FontSize', 14);
        subplot(2,2,2); comp2plot = dataFFT(:,dsearchn(Hz', sFreq)); topoplot(comp2plot./max(comp2plot),EEG.chanlocs,'maplimits',[-1 1],'numcontour',0,'electrodes','on','emarker2',{find(strcmpi({EEG.chanlocs.labels},electrode)) 'o' 'w' 4},'shading','interp');
        title([ 'Electrode Power at ' num2str(sFreq) ' Hz' ], 'FontSize', 14);
        subplot(2,2,[3:4]); plot(Hz,compSNR,'ro-','linew',1,'markersize',5,'markerface','w'); hold on;
        plot(Hz,elecSNR,'ko-','linew',1,'markersize',5,'markerface','w');
        set(gca,'xlim',xlim); xlabel('Frequency (Hz)', 'FontSize', 14), ylabel('SNR', 'FontSize', 14); legend({'Component'; electrode}, 'FontSize', 14); clear xlim
        saveas(figure(6), [directoryResults '/fig_ressVSelectrode.png']);

        %% Compute stability index

        % Gaussian filering extracted component
        compTimeConcat = reshape(compTime, 1,[]);
        compFiltered   = filterFGx(compTimeConcat, EEG.srate, sFreq, sFWHM);

        % Compute Hilbert Transform
        compHilbert = hilbert(compFiltered);

        % Extract phase angles
        compPhaseAngle = unwrap(angle(compHilbert));
        compPhaseHz    = (EEG.srate*diff(compPhaseAngle)) / (2*pi);

        % Apply a sliding moving median with a window width of 400 samples
        nOrder = 10;
        orders = linspace(10,400,nOrder)/2;
        orders = round(orders/(1000/EEG.srate));
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
        time = linspace(1, round(length(phaseMedFilt)/EEG.srate), length(phaseMedFilt));
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
        resultsEEG.(participantStr).([speedStr{iSpeed}]).stabilityIndex = stabilityIndex;
        resultsEEG.(participantStr).([speedStr{iSpeed}]).phaseAngle     = compPhaseAngle;
        resultsEEG.(participantStr).([speedStr{iSpeed}]).phaseHz        = phaseMedFilt;
        resultsEEG.(participantStr).([speedStr{iSpeed}]).compTime       = compTime;
        resultsEEG.(participantStr).([speedStr{iSpeed}]).Map            = compMaps(:,compMax);
        resultsEEG.(participantStr).([speedStr{iSpeed}]).chanLocs       = EEG.chanlocs;
        resultsEEG.(participantStr).([speedStr{iSpeed}]).compSNR        = compSNR;
        resultsEEG.(participantStr).([speedStr{iSpeed}]).elecSNR        = elecSNR;
        resultsEEG.(participantStr).([speedStr{iSpeed}]).compMax        = comp2Keep;
        resultsEEG.(participantStr).([speedStr{iSpeed}]).freqEEG        = EEG.srate;
        resultsEEG.(participantStr).([speedStr{iSpeed}]).freqStep       = sFreq;
        save([pathResults 'subAll/resultsEEG'], 'resultsEEG');

        clear comp2plot compFFT compFiltered compHilbert compMaps compPhaseAngle...
            compPhaseHz compSNR compTime compTimeConcat...
            dataTemp dataTempLength dataFFT elecSNR elecFFT Index...
            highCovariance highData Hz lambdaIndex Lambdas lambdaSorted lowCovariance lowData...
            phaseMed phaseMedFilt phaseTemp rCovariance sCovariance sData time timeVector W
        close all

    end

    clear dataAll Freq EEG elecLoc

end