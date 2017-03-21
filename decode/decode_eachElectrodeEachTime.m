function decode_eachElectrodeEachTime(data, savefilename, trialLabel,decodeSPGDir,times2downsample,freq,DP,C,nCond)
lastwarn('');
warning('on','all')

DP = decode_extractTimeFreqTrialInfo(DP,data);

clear P PPerm
isNotEnoughTrials = 0;

if DP.permute
    trialLabelPerm = zeros(DP.nPermutes,length(trialLabel));
    for iPerm = 1:DP.nPermutes
        trialLabelPerm(iPerm,:) = shuffle(trialLabel); % define shuffled trials before the rest of script. This keeps the time series intact
    end
end


for iTime = 1:DP.nTime
    if isNotEnoughTrials
        return
    end
    [tmp4 DP] = decode_extractFeatures(DP,data,iTime,[]);
    
    %% decode
    cResp = cell(nCond,1);
    for iCategory = 1:nCond
        cResp{iCategory} = [squeeze(tmp4(trialLabel == iCategory,:))];
        nTrialPerClass(iCategory) = size(cResp{iCategory},1);
    end

    if min(nTrialPerClass) < 2 ;
        disp(['# of trials are less than 2 for some class. skip.'])
        isNotEnoughTrials = 1;
        continue
    end
   
    [P{iTime}] = decodeRLSC(cResp,C,DP);

    [warnmsg, msgid] = lastwarn;
    if strcmp(msgid,'MATLAB:illConditionedMatrix')
        DP.warnmsg  = warnmsg;
        DP.msgid    = msgid;
        warning('off','last')
    end

    if ~mod(iTime,5)
        disp([ ': time = ' num2str(iTime) '/' num2str(DP.nTime) ' : ' datestr(now) ...
            ': # of trials = ' num2str(nTrialPerClass)])
    end
    
    %% permute decoding
    if DP.permute
        for iPerm = 1:DP.nPermutes
            cRespPerm = cell(nCond,1);
            for iCategory = 1:nCond
                switch DP.saveCategory
                    case 'pac'
                        % just shuffle power values
                        cRespPerm{iCategory}(:,1:DP.nFeaturesPhase)                 = [squeeze(tmp4(trialLabel == iCategory,1:DP.nFeaturesPhase))];                        
                        cRespPerm{iCategory}(:,DP.nFeaturesPhase+1:DP.nFeatures)    = [squeeze(tmp4(trialLabelPerm(iPerm,:) == iCategory,DP.nFeaturesPhase+1:DP.nFeatures))];                        
                    otherwise
                        cRespPerm{iCategory} = [squeeze(tmp4(trialLabelPerm(iPerm,:) == iCategory,:))];
                end
            end
            
            [Pperm{iTime,iPerm}] = decodeRLSC(cRespPerm,C,DP);
            
            if ~mod(iTime,5) && ~mod(iPerm,20)
                disp([ 'permute ' num2str(iPerm) ' of ' num2str(DP.nPermutes) ': time = ' num2str(iTime) '/' num2str(DP.nTime) ' : ' datestr(now) ...
                    ': # of trials = ' num2str(nTrialPerClass)])
            end
        end
    end
end


%% save for each lambda
% for iLam = 1:C.nLam

clear R Rperm
for iTime = 1:DP.nTime
    
    for iLambda = 1:C.nLam
        R.zmcorrTrain(iTime,:,iLambda)    = [P{iTime}(iLambda,:).zmcorrTrain];
        R.zmcorrTest(iTime,:,iLambda)     = [P{iTime}(iLambda,:).zmcorrTest];
        
        if DP.permute
            for iPerm = 1:DP.nPermutes
                Rperm.corrTrain(iTime,:,iLambda,iPerm) = [Pperm{iTime,iPerm}(iLambda,:).zmcorrTrain];
                Rperm.corrTest(iTime,:,iLambda,iPerm) = [Pperm{iTime,iPerm}(iLambda,:).zmcorrTest];
            end
        end
    end
end

if DP.permute
    save([decodeSPGDir 'dTF_' savefilename...
        '_' DP.ext '.mat'],...
        'R','Rperm','C','times2downsample','freq','nTrialPerClass','DP')
else
    save([decodeSPGDir 'dTF_' savefilename...
        '_' DP.ext '.mat'],...
        'R','C','times2downsample','freq','nTrialPerClass','DP')
end

if isdir('/gpfs/M2Home/projects/Monash052/jochem/')
    copyToKani(decodeSPGDir,['dTF_' savefilename '_' DP.ext '.mat'])
end

% end