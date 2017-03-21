function decodeAll(DEC,data,trialLabel,savefilename,decodeDir,times2downsample,freq)

%check if file already exists.
if exist([decodeDir 'dTF_' savefilename '.mat' ],'file'), disp(savefilename), return, end
disp(['working on: ' 'dTF_' savefilename '.mat'])

% get informaion about number of time and freq points etc.
DEC = decode_extractTimeFreqTrialInfo(DEC,data);

% extra check if matrix is composed the right way
if DEC.nTrials~=length(trialLabel)
    error('trial numbers don''t match')
end

if DEC.permute
    trialLabelPerm = zeros(DEC.nPermutes,length(trialLabel));
    for iPerm = 1:DEC.nPermutes
        trialLabelPerm(iPerm,:) = trialLabel(randperm(length(trialLabel))); % define shuffled trials before the rest of script. This keeps the time series intact
    end
else
    trialLabelPerm=[];
end

nTrialPerClass = NaN(1,DEC.nCond);
for iCategory = 1:DEC.nCond
    nTrialPerClass(iCategory) = length(find(trialLabel==iCategory));
end

decodePoint1 = [];
decodePoint2 = [];
switch DEC.dimension
    case {'eachElectrodeEachTime','multiElectrodeEachTime'};
        decodePoint1 = [DEC.nTime];
        decodePoint2 = 1;
    case {'eachElectrodeEachFreq','multiElectrodeEachFreq'};
        decodePoint1 = [DEC.nFreq];
        decodePoint2 = 1;
    case {'eachElectrodeEachTimeEachFreq','multiElectrodeEachTimeEachFreq'}
        decodePoint1 = [DEC.nTime];
        decodePoint2 = [DEC.nFreq];
    case {'eachElectrodeMultiTimeMultiFreq','multiElectrodeMultiTimeMultiFreq'}
        decodePoint1 = 1;
        decodePoint2 = 1;
end


% check whether we want to select a random subset of channels (in multi
% channel decoding). if so, loop over decoding script for number of
% iterations
if strfind(DEC.dimension, 'multi') && DEC.randChanSelection
    disp(['selecting ' num2str(DEC.randChanPercentage) '% of channels, with ' num2str(DEC.randChanCrossVal) ' iterations'])
    DEC.loop = DEC.randChanCrossVal;
else
    DEC.loop = 1;
end

% initialize variables
clear P PPerm
P = cell(DEC.loop,decodePoint1,decodePoint2);
Pperm = cell(DEC.loop,decodePoint1,decodePoint2);
decodeWarning = cell(DEC.loop,decodePoint1,decodePoint2);

isNotEnoughTrials = zeros(DEC.loop,decodePoint1,decodePoint2);
nFeat = zeros(DEC.loop,decodePoint1,decodePoint2);

% loop across nChannel subselections
for iLoop = 1:DEC.loop
    
    % select the random subset of channels on each iteration
    if DEC.randChanSelection
        DEC.nChan2select = ceil(DEC.randChanPercentage/100*DEC.nChan);
        DEC.chan2use = sort(randsample(DEC.nChan,DEC.nChan2select))';
        disp(['random channel iteration ' num2str(iLoop) ' of ' num2str(DEC.randChanCrossVal)])
    else
        DEC.chan2use = 1:DEC.nChan;
    end
    
    %% loop across time and frequency points
    for iDecPoints1 = 1:decodePoint1
        for iDecPoints2 = 1:decodePoint2
            
            % decode the time and frequency point
            [P{iLoop,iDecPoints1,iDecPoints2}, Pperm{iLoop,iDecPoints1,iDecPoints2}, decodeWarning{iLoop,iDecPoints1,iDecPoints2}, isNotEnoughTrials(iLoop,iDecPoints1,iDecPoints2), nFeat(iLoop,iDecPoints1,iDecPoints2)] = decodeDecPoint(DEC,data,trialLabel,trialLabelPerm,iDecPoints1,iDecPoints2);
            
        end
        if ~mod(iDecPoints1,5)
            disp([ ': loop1 = ' num2str(iDecPoints1) '/' num2str(decodePoint1) ' : ' datestr(now) ': # of trials = ' num2str(nTrialPerClass)])
        end
    end
end

DEC.nFeat = nFeat;
if ~isempty(find(isNotEnoughTrials,1))
    return
end
decodeWarning;

%% extract decode accuracies

% extract decoding accuracies from structures/cells
clear R Rperm
for iLoop = 1:DEC.loop
    
    for iDecPoints1 = 1:decodePoint1
        for iDecPoints2 = 1:decodePoint2
            for iLambda = 1:DEC.nLam
                tmpR.zmcorrTrain(iLoop,iDecPoints1,iDecPoints2,:,iLambda)    = [P{iLoop,iDecPoints1,iDecPoints2}(iLambda,:).zmcorrTrain];
                tmpR.zmcorrTest(iLoop,iDecPoints1,iDecPoints2,:,iLambda)     = [P{iLoop,iDecPoints1,iDecPoints2}(iLambda,:).zmcorrTest];
                
                
                if DEC.permute
                    for iPerm = 1:DEC.nPermutes
                        tmpRperm.zmcorrTrain(iLoop,iDecPoints1,iDecPoints2,:,iLambda,iPerm) = [Pperm{iLoop,iDecPoints1,iDecPoints2}{iPerm}(iLambda,:).zmcorrTrain];
                        tmpRperm.zmcorrTest(iLoop,iDecPoints1,iDecPoints2,:,iLambda,iPerm) = [Pperm{iLoop,iDecPoints1,iDecPoints2}{iPerm}(iLambda,:).zmcorrTest];
                        
                    end
                end
            end
        end
    end
end

% mean across random channel subselections, dimension reduction 
for iDecPoints1 = 1:decodePoint1
    for iDecPoints2 = 1:decodePoint2
        for iLambda = 1:DEC.nLam
            R.zmcorrTrain(iDecPoints1,iDecPoints2,:,iLambda)       = mean(tmpR.zmcorrTrain(:,iDecPoints1,iDecPoints2,:,iLambda),1);
            R.zmcorrTest(iDecPoints1,iDecPoints2,:,iLambda)        = mean(tmpR.zmcorrTest(:,iDecPoints1,iDecPoints2,:,iLambda),1);
            R.stdZmcorrTrain(iDecPoints1,iDecPoints2,:,iLambda)    = std(tmpR.zmcorrTrain(:,iDecPoints1,iDecPoints2,:,iLambda),1);
            R.stdZmcorrTest(iDecPoints1,iDecPoints2,:,iLambda)     = std(tmpR.zmcorrTest(:,iDecPoints1,iDecPoints2,:,iLambda),1);
        end
    end
end

if DEC.permute
    for iDecPoints1 = 1:decodePoint1
        for iDecPoints2 = 1:decodePoint2
            for iLambda = 1:DEC.nLam
                Rperm.zmcorrTrain(iDecPoints1,iDecPoints2,:,iLambda,:)   = mean(tmpRperm.zmcorrTrain(:,iDecPoints1,iDecPoints2,:,iLambda,:),1);
                Rperm.zmcorrTest(iDecPoints1,iDecPoints2,:,iLambda,:)    = mean(tmpRperm.zmcorrTest(:,iDecPoints1,iDecPoints2,:,iLambda,:),1);
            end
        end
    end
end

% save results
try
    if DEC.permute
        save([decodeDir 'dTF_' savefilename '.mat'],...
            'R','Rperm','times2downsample','freq','nTrialPerClass','DEC','decodeWarning')
    else
        save([decodeDir 'dTF_' savefilename '.mat'],...
            'R','times2downsample','freq','nTrialPerClass','DEC','decodeWarning')
    end
catch
    disp(['could not save ' savefilename])
    keyboard
end

if isdir('/gpfs/M2Home/projects/Monash052/jochem/')
    copyToKani(decodeSPGDir,['dTF_' savefilename '.mat'])
end

if isdir('/gpfs/M1Scratch/Monash052/')
    copyToKani_M2(decodeSPGDir,['dTF_' savefilename '.mat'])
end


end

%% extra functions


function [P, Pperm, decodeWarning, isNotEnoughTrials, nFeat] = decodeDecPoint(DEC,data,trialLabel,trialLabelPerm,iDecPoint1,iDecPoint2)
tmpData=[];

% extract relevant features
switch DEC.dimension
    case {'eachElectrodeEachTime','multiElectrodeEachTime'};
        [tmpData, nFeat] = decode_extractFeatures(DEC,data,iDecPoint1);% only time is a decode point
    case {'eachElectrodeEachFreq','multiElectrodeEachFreq'};
        [tmpData, nFeat] = decode_extractFeatures(DEC,data,[],iDecPoint1);% only freq is a decode point
    case {'eachElectrodeEachTimeEachFreq','multiElectrodeEachTimeEachFreq'}
        [tmpData, nFeat] = decode_extractFeatures(DEC,data,iDecPoint1,iDecPoint2);% time and freq are decode points
    case {'eachElectrodeMultiTimeMultiFreq','multiElectrodeMultiTimeMultiFreq'}
        [tmpData, nFeat] = decode_extractFeatures(DEC,data);% time and freq are used as separate features in one big 'decode point'
end

% decode
if DEC.permute
    [P, Pperm,  decodeWarning, isNotEnoughTrials] = decodeFeature_RLSC(tmpData,DEC,trialLabel,trialLabelPerm,DEC.nCond);
else
    [P, ~,      decodeWarning, isNotEnoughTrials] = decodeFeature_RLSC(tmpData,DEC,trialLabel,[],DEC.nCond);
    Pperm = cell(1);
end
end

