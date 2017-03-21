function [signLocalizerChannels localizerChannelsDim decodeChan times2downsample] = getSignificantIowaLocalizerChannels(subID,SPG,DP,decodeCheckDir,EXP,threshold)
subSpecsIowaFaceLocalizer
% setDirIowaFaceLocalizer

% tmpLocalizerChannels = dir([subDir filesep 'checkDecodePerformance' filesep '*'  Job.localizerChannelsCondition '*']);
% 
% % filenameLocalizer = ['dTF_' subID '_li' num2str(1) '_distrib_4_ses_' Job.localizerChannelsCondition '_' DP.ext '.mat'];
% localizerChannels = zeros(1,length(tmpLocalizerChannels));
% for iChan = 1:length(tmpLocalizerChannels)
%     localizerChannels(iChan)  = str2double(tmpLocalizerChannels(iChan).name(length(['dTF_' subID '_li'])+1:...
%         (length(tmpLocalizerChannels(iChan).name)-length(['_distrib_4_ses_' Job.localizerChannelsCondition '_' DP.ext '.mat']))));
% end
% localizerChannels = sort(localizerChannels);

% decodeChan = struct([]);

if isempty(SPG) || SPG.nSessionDecode==NaN
    SPG.nSessionDecode=1;
end
if isempty(DP)
    DP.vDim = [1 2];
end
if isempty(decodeCheckDir)

end


signLocalizerChannels    = zeros(size(Channels2check));
localizerChannelsDim     = zeros(length(Channels2check),4);
for iChan = Channels2check
    val = ['li' num2str(iChan)];
    filename = ['dTF_' subID '_' val '_distrib_' num2str(SPG.nSessionDecode) '_ses_' EXP.cond '_' DP.ext];

    if ~exist([decodeCheckDir filesep filename '.mat'],'file')
        filename
        keyboard
    end
    decodeChan(iChan) = load([decodeCheckDir filesep filename '.mat']);
    
    for iDim = DP.vDim
        
        
        if isempty(threshold)
            if ~isempty(decodeChan(iChan).pCorr(iDim).pCorr)
                signLocalizerChannels(iChan) = 1;
                localizerChannelsDim(iChan,iDim) = 1;
            end
        else 
            if ~isempty(find(squeeze(decodeChan(iChan).pMeanTestCorr(iDim,:,:))>threshold))
                signLocalizerChannels(iChan) = 1;
                localizerChannelsDim(iChan,iDim) = 1;
            end
            
        end
    end
end

times2downsample = decodeChan(1).times2downsample;