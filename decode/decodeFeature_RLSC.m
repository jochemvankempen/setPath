function [P, Pperm, decodeWarning, isNotEnoughTrials] = decodeFeature_RLSC(tmp4,DEC,trialLabel,trialLabelPerm,nCond)


lastwarn('');
warning('on','all')

%% decode
cResp = cell(nCond,1);
for iCategory = 1:nCond
    cResp{iCategory} = [squeeze(tmp4(trialLabel == iCategory,:))];
    nTrialPerClass(iCategory) = size(cResp{iCategory},1);
end
isNotEnoughTrials=0;
if min(nTrialPerClass) < 4 ;
    disp(['# of trials are less than 3 for some class. skip.'])
    P=[];
    Pperm=[];
    decodeWarning=[];
    isNotEnoughTrials=1;
    return
end

[P] = decodeRLSC(cResp,DEC);
[warnmsg, msgid] = lastwarn;
if strcmp(msgid,'MATLAB:illConditionedMatrix')
    decodeWarning.warnmsg  = warnmsg;
    decodeWarning.msgid    = msgid;
%     keyboard
%     warning('off','last')
else
    decodeWarning.warnmsg  = num2str(NaN);
    decodeWarning.msgid    = num2str(NaN);
end

%% permute decoding
if DEC.permute
    parfor iPerm = 1:DEC.nPermutes
%         clear tmpPerm
        cRespPerm = cell(nCond,1);
        %                 for iCategory = 1:nCond
        %                     cRespPerm{iCategory} = [squeeze(tmp4(vPosCorrTrialsPerm(iPerm,:) == iCategory,:))];
        %                 end
        for iCategory = 1:nCond
            switch DEC.saveCategory
                case 'pac'
                    % just shuffle power values
                    cRespPerm{iCategory}(:,1:DEC.nFeaturesPhase)                 = [squeeze(tmp4(trialLabel == iCategory,1:DEC.nFeaturesPhase))];
                    cRespPerm{iCategory}(:,DEC.nFeaturesPhase+1:DEC.nFeatures)    = [squeeze(tmp4(trialLabelPerm(iPerm,:) == iCategory,DEC.nFeaturesPhase+1:DEC.nFeatures))];
                otherwise
                    cRespPerm{iCategory} = [squeeze(tmp4(trialLabelPerm(iPerm,:) == iCategory,:))];
            end
        end
        
        Pperm{iPerm} = decodeRLSC(cRespPerm,DEC);
        if ~mod(iPerm,20)
%             
            disp([ 'permute ' num2str(iPerm) ' of ' num2str(DEC.nPermutes) ' : ' datestr(now) ...
                ': # of trials = ' num2str(nTrialPerClass)])
        end
        
    end
else
    Pperm=[];
end
Pperm = Pperm';
end