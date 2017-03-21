% load(['/Volumes/backupCFSBMDF/IowaData/MAP/hhbcntInfo.mat'])
% load([mapDir '/bcntInfo.mat'])
% load([mapDir '/bcntInfo.mat'])


% CNT = [];
switch subID
    case '146'
        cntID = '146_24_CFS';
    case '147'
        cntID = '147_DF_CFS';
    case '154a'
        cntID = '154_CFS_73_74';
    case '154b'
        cntID = '154_DF_BM_CFS_137';
    case '154'
        cntID = '154_combined';
    case '156'
        cntID = '156_CFS';
    otherwise
        cntID = subID;
end


try
    load([DIR.map 'bcntInfo_' cntID '.mat'],'bCNT')
catch
    disp('error : bcntInfo not found')
    return
end

CNT=bCNT;

% for iTmp = 1:length(bcntInfo)
%     if strcmp( bcntInfo(iTmp).id , cntID)
%         CNT = bcntInfo(iTmp);
%         break
%     end
% end
% if isempty(CNT)
%     disp('error : bcntInfo not found')
%     keyboard
% end
%%
nElectrode = length(CNT.label);
disp(['# of bipolar electrodes = ' num2str(nElectrode)])
if nElectrode < 170
    disp('too few electrodes')
    keyboard
end
%% check bipolar locations
% 040313 bipolar electrode locations are now checked directly after setting vElectrode in phi_Comp_SegmentedData*.m
check = 0;
if check
    tmpBrainView = {'temporal','ventral'};
    getContactInfoForMappping
    vBrainView = 1:2;
    disp('check bipolar electrode mapping in loadContactInfo_bipolar_CFS')
    for iBrainView = vBrainView
        alignedMap{iBrainView} = load([mapMatFileDir '/' subID '/aligned_' tmpBrainView{iBrainView} '.mat']);
        figure(iBrainView),clf
        imagesc(uint8(alignedMap{iBrainView}.img))
        hold on
        axis image
        
        %%
        clear pairBipolar
        for iBipolar = cBipolar{iBrainView}
            iBipolar
            cBipolar{iBrainView}(1)
            jBipolar = iBipolar-cBipolar{iBrainView}(1)+1;
            tmp = CNT.label{iBipolar};
            pairBipolar(jBipolar,1) = str2num(tmp(end-6:end-4));
            pairBipolar(jBipolar,2) = str2num(tmp(end-2:end));
            disp([CNT.label{iBipolar} ' : pairBipolar =' num2str(pairBipolar(jBipolar,:))])
            
            pairBipolar2(jBipolar,:) = pairBipolar(jBipolar,:) - cElectrodes{iBrainView}(1) + 1;
            
            tmpx = alignedMap{iBrainView}.xcent(pairBipolar2(jBipolar,:));
            tmpy = alignedMap{iBrainView}.ycent(pairBipolar2(jBipolar,:));
            plot(tmpx,tmpy,'b-','linewidth',3)
        end
        print(gcf,'-dpng',[mapMatFileDir '/' subID '/checkBipolar_' tmpBrainView{iBrainView} '.png'])
        
    end
end

%%
check2 = 0;
if check2
    %% check electrodes
    figure(100),clf
    imagesc(uint8(alignedMap{iBrainView}.img))
    axis image
    hold on
    %
    for iElectrode = cCheckKeyElectrodes{iBrainView}
        %% map onto the map
        jElectrode = iElectrode-cElectrodes{iBrainView}(1)+1;
        x = alignedMap{iBrainView}.xcent(jElectrode);
        y = alignedMap{iBrainView}.ycent(jElectrode);
        plot(x,y,'o')
        fs = 10
        text(x,y,CNT.label{iElectrode},'interpret','none','fontsize',fs)
    end
    %print(gcf,'-dpng',[mapMatFileDir '/' subID '/checkKeyElectrode_' tmpBrainView{iBrainView} '.png'])
end