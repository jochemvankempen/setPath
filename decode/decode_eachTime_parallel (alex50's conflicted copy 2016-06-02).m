function decode_eachTime_parallel(data, savefilename, trialLabel,decodeSPGDir,times2downsample,freq,DEC,nCond)

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
    trialLabelPerm=[];
end


P = cell(1,DEC.nTime);
Pperm = cell(1,DEC.nTime);
%Pperm = cell(DEC.nTime,DEC.nPermutes);
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
    
    parfor iTime = 1:DEC.nTime
        tmp4=[];
%         decodeWarning=[];
        %     if isNotEnoughTrials(iTime)
        %         return
        %     end
        [tmp4, nFeat{iLoop,iTime}] = decode_extractFeatures(DEC,data,iTime);
        
        if DEC.permute
            [P{iLoop,iTime}, Pperm{iLoop,iTime}, decodeWarning{iLoop,iTime}, isNotEnoughTrials(iLoop,iTime)] = decodeFeature_RLSC(tmp4,DEC,trialLabel,trialLabelPerm,nCond);
        else
            [P{iLoop,iTime}, ~, decodeWarning{iLoop,iTime}, isNotEnoughTrials(iLoop,iTime)] = decodeFeature_RLSC(tmp4,DEC,trialLabel,[],nCond);
        end
        %%
        
        if ~mod(iTime,5)
            disp([ ': time = ' num2str(iTime) '/' num2str(DEC.nTime) ' : ' datestr(now) ...
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
    
    for iTime = 1:DEC.nTime
        
        for iLambda = 1:DEC.nLam
            tmpR.zmcorrTrain(iLoop,iTime,:,iLambda)    = [P{iLoop,iTime}(iLambda,:).zmcorrTrain];
            tmpR.zmcorrTest(iLoop,iTime,:,iLambda)     = [P{iLoop,iTime}(iLambda,:).zmcorrTest];
            
            
            if DEC.permute
                for iPerm = 1:DEC.nPermutes
                    tmpRperm.zmcorrTrain(iLoop,iTime,:,iLambda,iPerm) = [Pperm{iLoop,iTime}{iPerm}(iLambda,:).zmcorrTrain];
                    tmpRperm.zmcorrTest(iLoop,iTime,:,iLambda,iPerm) = [Pperm{iLoop,iTime}{iPerm}(iLambda,:).zmcorrTest];

                end
            end
        end
    end
end

for iTime = 1:DEC.nTime
    for iLambda = 1:DEC.nLam
        R.zmcorrTrain(iTime,:,iLambda)       = mean(tmpR.zmcorrTrain(:,iTime,:,iLambda),1);
        R.zmcorrTest(iTime,:,iLambda)        = mean(tmpR.zmcorrTest(:,iTime,:,iLambda),1);
        R.stdZmcorrTrain(iTime,:,iLambda)    = std(tmpR.zmcorrTrain(:,iTime,:,iLambda),1);
        R.stdZmcorrTest(iTime,:,iLambda)     = std(tmpR.zmcorrTest(:,iTime,:,iLambda),1);
    end
end
        
% R.zmcorrTrain       = squeeze(mean(tmpR.zmcorrTrain,1));
% R.zmcorrTest        = squeeze(mean(tmpR.zmcorrTest,1));
% R.stdZmcorrTrain    = squeeze(std(tmpR.zmcorrTrain,1));
% R.stdZmcorrTest     = squeeze(std(tmpR.zmcorrTest,1));
              
if DEC.permute
    for iTime = 1:DEC.nTime
        for iLambda = 1:DEC.nLam
            Rperm.zmcorrTrain(iTime,:,iLambda,:)   = mean(tmpRperm.zmcorrTrain(:,iTime,:,iLambda,:),1);
            Rperm.zmcorrTest(iTime,:,iLambda,:)    = mean(tmpRperm.zmcorrTest(:,iTime,:,iLambda,:),1);
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