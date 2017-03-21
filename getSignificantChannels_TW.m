function [signChannels signDim signMW decodeChan times2downsample signChannelsAllDimAllTW thresholdSetting] = getSignificantChannels_TW(subID,arrayID,SPG,DP,EXP)
subSpecsIowaFaceLocalizer
% setDirIowaFaceLocalizer
%%

% thresholdSetting = 'allTimepoints';
% thresholdSetting = 'threshold';
thresholdSetting = 'consecutiveTimepoints';

switch thresholdSetting
    case 'allTimepoints'
        disp('using all significant channels, single timepoint, one dimension')
    case 'threshold'
        threshold = 0.3;
        disp(['using threshold ' num2str(threshold)])
    case 'consecutiveTimepoints'
        consecutiveTimepoints = 3;
        disp(['using threshold, multiple (n=' num2str(consecutiveTimepoints) ') consecutive timepoints'])
        
end


%%
signChannels    = zeros(1,Channels2check(end));
signChannelsAllDimAllTW    = zeros(1,Channels2check(end));
signDim     = zeros(Channels2check(end),length(DP.vDim));
signMW      = zeros(Channels2check(end),length(DP.vDim),length(SPG.TW),length(SPG.TW));
for iChan = (Channels2check(length(Channels2checkOriginal)+1:end))
    val = ['li' num2str(iChan)];
    filename = ['dTF_' subID '_' val '_movingwin_' EXP.cond '_' StimCategory '_' DP.decodeSPG '_' DP.ext];
    
    if ~exist([DP.decodeCheckDir filesep filename '.mat'],'file')
        filename
        keyboard
        continue
    end
    decodeChan(iChan) = load([DP.decodeCheckDir filesep filename '.mat']);
end
%%



for iChan = (Channels2check(length(Channels2checkOriginal)+1:end))
    switch DP.dimension
        case 'eachElectrodeEachTime'
            for iDim = DP.vDim
                for iTW =DP.nMovingwin
                    %                     if ~DP.movingWindowSizes2use(iTW)
                    %                         continue
                    %                     end
                    
                    
                    for iTTW =DP.nTotalTimewin
                        if iTTW < iTW
                            continue
                        end
                        switch thresholdSetting
                            case 'allTimepoints'
                                if ~isempty(decodeChan(iChan).pCorr(iDim,iTW,iTTW).pCorr)
                                    signChannels    (iChan)                                 = 1;
                                    signDim         (iChan,iDim)                            = 1;
                                    signMW          (iChan,iDim,iTW,iTTW)   = 1;
                                    
                                    %                                 signAllMW
                                end
                            case 'threshold'
                                
                                if ~isempty(find(decodeChan(iChan).pMeanTestCorr(iDim,iTW,iTTW,:)>threshold))
                                    signChannels    (iChan)                                 = 1;
                                    signDim         (iChan,iDim)                            = 1;
                                    signMW          (iChan,iDim,iTW,iTTW)   = 1;
                                end
                            case 'consecutiveTimepoints'
                                
                                if length(decodeChan(iChan).pCorr(iDim,iTW,iTTW).pCorr)<2
                                    continue
                                else
                                    switch consecutiveTimepoints
                                        case 2
                                            for iSignPoint = 2:length(decodeChan(iChan).pCorr(iDim,iTW,iTTW).pCorr)
                                                if decodeChan(iChan).pCorr(iDim,iTW,iTTW).pCorr(iSignPoint) - decodeChan(iChan).pCorr(iDim,iTW,iTTW).pCorr(iSignPoint-1) == 1
                                                    
                                                    signChannels    (iChan)                                 = 1;
                                                    signDim         (iChan,iDim)                            = 1;
                                                    signMW          (iChan,iDim,iTW,iTTW)   = 1;
                                                end
                                            end
                                        case 3
                                            for iSignPoint = 3:length(decodeChan(iChan).pCorr(iDim,iTW,iTTW).pCorr)
                                                if decodeChan(iChan).pCorr(iDim,iTW,iTTW).pCorr(iSignPoint) - decodeChan(iChan).pCorr(iDim,iTW,iTTW).pCorr(iSignPoint-1) == 1 && ...
                                                        decodeChan(iChan).pCorr(iDim,iTW,iTTW).pCorr(iSignPoint-1) - decodeChan(iChan).pCorr(iDim,iTW,iTTW).pCorr(iSignPoint-2) == 1
                                                    
                                                    signChannels    (iChan)                                 = 1;
                                                    signDim         (iChan,iDim)                            = 1;
                                                    signMW          (iChan,iDim,iTW,iTTW)   = 1;
                                                end
                                            end
                                        otherwise
                                            keyboard
                                    end
                                end
                        end
                    end
                    
                end
            end
            
            
            %             for iTW = DP.nMovingwin
            %
            %
            %                 for iTTW = DP.nTotalTimewin
            if length(find(squeeze(signMW(iChan,:,DP.nMovingwin,DP.nTotalTimewin))>0.3)) == 20
                signChannelsAllDimAllTW(iChan) =1;
            end
            %                 end
            %             end
        case 'eachElectrodeEachTimeEachFreq'
            for iDim = DP.vDim
                for iTW =DP.nMovingwin
                    for iTTW =DP.nTotalTimewin
                        if iTTW < iTW
                            continue
                        end
                        if 0
                            if ~isempty(decodeChan(iChan).pCorr(iDim,iTW,iTTW).pCorr)
                                signChannels    (iChan)                                 = 1;
                                signDim         (iChan,iDim)                            = 1;
                                signMW          (iChan,iDim,iTW,iTTW)   = 1;
                            end
                        end
                        if ~isempty(find(decodeChan(iChan).pMeanTestCorr(iDim,iTW,iTTW,:,:)>0.3))
                            signChannels    (iChan)                                 = 1;
                            signDim         (iChan,iDim)                            = 1;
                            signMW          (iChan,iDim,iTW,iTTW)   = 1;
                        end
                    end
                end
            end
            
    end
end

times2downsample = decodeChan(arrayID(1)).times2downsample;