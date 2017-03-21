function decode_eachFreq_parallel(data, savefilename, trialLabel,decodeSPGDir,times2downsample,freq,DEC,nCond)

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
else
    trialLabelPerm = [];
end


P = cell(1,DEC.nFreq);
Pperm = cell(1,DEC.nFreq);
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
    
    parfor iFreq = 1:DEC.nFreq
        tmp4=[];
        %     if isNotEnoughTrials(iFreq)
        %         return
        %     end

        [tmp4 nFeat{iLoop,iFreq}] = decode_extractFeatures(DEC,data,[],iFreq);
        
        if DEC.permute
            [P{iLoop,iFreq}, Pperm{iLoop,iFreq}, decodeWarning{iLoop,iFreq}, isNotEnoughTrials(iLoop,iFreq)] = decodeFeature_RLSC(tmp4,DEC,trialLabel,trialLabelPerm,nCond);
        else
            [P{iLoop,iFreq}, ~, decodeWarning{iLoop,iFreq}, isNotEnoughTrials(iLoop,iFreq)] = decodeFeature_RLSC(tmp4,DEC,trialLabel,[],nCond);
        end

        if ~mod(iFreq,5)
            disp([ ': freq = ' num2str(iFreq) '/' num2str(DEC.nFreq) ' : ' datestr(now) ...
                ': # of trials = ' num2str(nTrialPerClass)])
        end
    end
end
DEC.nFeat = nFeat;
if ~isempty(find(isNotEnoughTrials))
    return
end
decodeWarning;
%% save for each lambda
% for iLam = 1:DEC.nLam

clear R Rperm
for iLoop = 1:DEC.loop
    
    for iFreq = 1:DEC.nFreq
        
        for iLambda = 1:DEC.nLam
            tmpR.zmcorrTrain(iLoop,iFreq,:,iLambda)    = [P{iLoop,iFreq}(iLambda,:).zmcorrTrain];
            tmpR.zmcorrTest(iLoop,iFreq,:,iLambda)     = [P{iLoop,iFreq}(iLambda,:).zmcorrTest];
            
            
            if DEC.permute
                for iPerm = 1:DEC.nPermutes
                    tmpRperm.zmcorrTrain(iLoop,iFreq,:,iLambda,iPerm) = [Pperm{iLoop,iFreq}{iPerm}(iLambda,:).zmcorrTrain];
                    tmpRperm.zmcorrTest(iLoop,iFreq,:,iLambda,iPerm) = [Pperm{iLoop,iFreq}{iPerm}(iLambda,:).zmcorrTest];

                end
            end
        end
    end
end

for iFreq = 1:DEC.nFreq
    for iLambda = 1:DEC.nLam
        R.zmcorrTrain(iFreq,:,iLambda)       = mean(tmpR.zmcorrTrain(:,iFreq,:,iLambda),1);
        R.zmcorrTest(iFreq,:,iLambda)        = mean(tmpR.zmcorrTest(:,iFreq,:,iLambda),1);
        R.stdZmcorrTrain(iFreq,:,iLambda)    = std(tmpR.zmcorrTrain(:,iFreq,:,iLambda),1);
        R.stdZmcorrTest(iFreq,:,iLambda)     = std(tmpR.zmcorrTest(:,iFreq,:,iLambda),1);
    end
end
        
% R.zmcorrTrain       = squeeze(mean(tmpR.zmcorrTrain,1));
% R.zmcorrTest        = squeeze(mean(tmpR.zmcorrTest,1));
% R.stdZmcorrTrain    = squeeze(std(tmpR.zmcorrTrain,1));
% R.stdZmcorrTest     = squeeze(std(tmpR.zmcorrTest,1));
              
if DEC.permute
    for iFreq = 1:DEC.nFreq
        for iLambda = 1:DEC.nLam
            Rperm.zmcorrTrain(iFreq,:,iLambda,:)   = mean(tmpRperm.zmcorrTrain(:,iFreq,:,iLambda,:),1);
            Rperm.zmcorrTest(iFreq,:,iLambda,:)    = mean(tmpRperm.zmcorrTest(:,iFreq,:,iLambda,:),1);
        end
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
catch
    disp(['could not save ' savefilename])
end
if isdir('/gpfs/M2Home/projects/Monash052/jochem/')
    copyToKani(decodeSPGDir,['dTF_' savefilename '.mat'])
end

if isdir('/gpfs/M1Scratch/Monash052/')
    copyToKani_M2(decodeSPGDir,['dTF_' savefilename '.mat'])
end
% end