function decode_multiLoopMultiFreq(data, savefilename, trialLabel,decodeSPGDir,times2downsample,freq,DEC,nCond)

DEC = decode_extractTimeFreqTrialInfo(DEC,data);
if DEC.nTrials~=length(trialLabel)
    error('trial numbers don''t match')
end
clear P PPerm
isNotEnoughTrials = 0;

if DEC.permute
    trialLabelPerm = zeros(DEC.nPermutes,length(trialLabel));
    for iPerm = 1:DEC.nPermutes
        trialLabelPerm(iPerm,:) = trialLabel(randperm(length(trialLabel))); % define shuffled trials before the rest of script. This keeps the time series intact
    end
end


P = cell(1,DEC.nTime);
Pperm = cell(1,DEC.nTime);
for iCategory = 1:nCond
    nTrialPerClass(iCategory) = length(find(trialLabel==iCategory));
end


if strfind(DEC.dimension, 'multi') && DEC.randChanSelection
    disp(['selecting ' num2str(DEC.randChanPercentage) '% of channels, with ' num2str(DEC.randChanCrossVal) ' iterations'])
    DEC.loop = DEC.randChanCrossVal;
else 
    DEC.loop = 1;
end


for iLoop = 1:DEC.loop
    if DEC.randChanSelection
        DEC.nChan2select = ceil(DEC.randChanPercentage/100*DEC.nChan);
        DEC.chan2use = sort(randsample(DEC.nChan,DEC.nChan2select))';
        disp(['random channel iteration ' num2str(iLoop) ' of ' num2str(DEC.randChanCrossVal)])
    else
        DEC.chan2use = 1:DEC.nChan;
    end
    
    
    [tmp4 nFeat(iLoop)] = decode_extractFeatures(DEC,data);
    
    if DEC.permute
        [P{iLoop}, Pperm{iLoop}, decodeWarning(iLoop).warn isNotEnoughTrials(iLoop)] = decodeFeature_RLSC(tmp4,DEC,trialLabel,trialLabelPerm,nCond);
    else
        [P{iLoop}, ~, decodeWarning(iLoop).warn isNotEnoughTrials(iLoop)] = decodeFeature_RLSC(tmp4,DEC,trialLabel,[],nCond);
    end
end
%%

DEC.nFeat = nFeat;
if ~isempty(find(isNotEnoughTrials))
    return
end

%% save for each lambda
% for iLam = 1:DEC.nLam

clear R Rperm
for iLoop = 1:DEC.loop
    
    
    for iLambda = 1:DEC.nLam
        tmpR.zmcorrTrain(iLoop,:,iLambda)    = [P{iLoop}(iLambda,:).zmcorrTrain];
        tmpR.zmcorrTest(iLoop,:,iLambda)     = [P{iLoop}(iLambda,:).zmcorrTest];
        
        if DEC.permute
            for iPerm = 1:DEC.nPermutes
                tmpRperm.zmcorrTrain(iLoop,:,iLambda,iPerm) = [Pperm{iLoop}{iPerm}(iLambda,:).zmcorrTrain];
                tmpRperm.zmcorrTest(iLoop,:,iLambda,iPerm) = [Pperm{iLoop}{iPerm}(iLambda,:).zmcorrTest];
            end
        end
    end
end

for iLambda = 1:DEC.nLam % to keep format the same as other decoding scripts
    R.zmcorrTrain(1,:,iLambda)       = mean(tmpR.zmcorrTrain(:,:,iLambda),1);
    R.zmcorrTest(1,:,iLambda)        = mean(tmpR.zmcorrTest(:,:,iLambda),1);
    R.stdZmcorrTrain(1,:,iLambda)    = std(tmpR.zmcorrTrain(:,:,iLambda),1);
    R.stdZmcorrTest(1,:,iLambda)     = std(tmpR.zmcorrTest(:,:,iLambda),1);
end
                      
if DEC.permute
    for iLambda = 1:DEC.nLam
        Rperm.zmcorrTrain(1,:,iLambda,:)   = mean(tmpRperm.zmcorrTrain(:,:,iLambda,:),1);
        Rperm.zmcorrTest(1,:,iLambda,:)    = mean(tmpRperm.zmcorrTest(:,:,iLambda,:),1);
    end
end



try
    if DEC.permute
        save([decodeSPGDir 'dTF_' savefilename '.mat'],...
            'R','Rperm','times2downsample','freq','nTrialPerClass','DEC','decodeWarning')
    else
        save([decodeSPGDir 'dTF_' savefilename '.mat'],...
            'R','times2downsample','freq','nTrialPerClass','DEC','decodeWarning')
    end
catch me
    disp(['could not save ' savefilename])
end

if isdir('/gpfs/M2Home/projects/Monash052/jochem/')
    copyToKani(decodeSPGDir,['dTF_' savefilename '.mat'])
end

if isdir('/gpfs/M1Scratch/Monash052/')
    copyToKani_M2(decodeSPGDir,['dTF_' savefilename '.mat'])
end
% end