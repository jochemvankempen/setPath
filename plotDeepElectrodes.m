function plotDeepElectrodes(subID,PLOT,savefilename,figDir,channelValues,colorScale)
% plots all labels none electrodes
subSpecsIowaFaceLocalizer

if isempty(colorScale)
    scaleDecode = [min(channelValues) max(channelValues)];
else
    scaleDecode = colorScale;
end

endName = 5;
getStripsInfo
spacing_yAxis = zeros(size(labels_none));
nChan2use    = zeros(size(labels_none));
%% 
%%
figureSettings
axis off
hold on
r = 0.03;
theta = linspace(0,2*pi,16);
R = r*[cos(theta); sin(theta)];


for iLabel = 1:length(labels_none)
    
    deepChanIdx     = [];
    deepBipChanIdx  = [];
    
    for iChan = 1:length(bCNT.label)
        
        if strcmpi(bCNT.label{iChan}(1:end-(endName)),labels_none(iLabel))
            deepChanIdx = [deepChanIdx iChan];
        elseif strcmpi([bCNT.label{iChan}(2:end-(endName*2))],labels_none(iLabel))
            deepBipChanIdx = [deepBipChanIdx iChan];
        end
    end

    switch PLOT.channels
        case 'unipolar'
            deepChan2use{iLabel} = deepChanIdx;
        case 'bipolar'
            deepChan2use{iLabel} = deepBipChanIdx;
        case 'all'
            deepChan2use{iLabel} = [deepChanIdx deepBipChanIdx];
    end
    
    nChan2use(iLabel) = length(deepChan2use{iLabel});
    
    spacing_yAxis(iLabel) = 0+iLabel*0.8/length(labels_none);

    text(0,spacing_yAxis(iLabel),labels_none{iLabel},'interpret','none','FontSize',15)

end

spacing_xAxis = .5;%5/max(nChan2use);

clear colormap cmap
colormap('default');

for iLabel = 1:length(labels_none)
    
    for iChan = 1:1:length(deepChan2use{iLabel})
%                 colval = (channelValues(deepChan2use{iLabel}(iChan))-scaleDecode(1))/diff(scaleDecode)*255;
        colval = (channelValues(deepChan2use{iLabel}(iChan)));
        patch(spacing_xAxis*iChan+(R(1,:)*5)+1,spacing_yAxis(iLabel)+R(2,:)-(1/length(labels_none)/4),colval)
        
    end   
    
end
xlim([0 spacing_xAxis*max(nChan2use)+1.5])
hold off
% h = colorbar('peer',gca);
% set(h,'ylim',[0 255],'ytick',linspace(0,255,5),'yticklabel',num2str(linspace(scaleDecode(1),scaleDecode(2),5)',2),'FontSize',15,'colormap','Jet');
% set(get(h,'xlabel'),'String', 'A''','FontSize',20);

%%
caxis([scaleDecode(1) scaleDecode(2)]);
cbh=colorbar;
set(cbh,'YTick',linspace(scaleDecode(1),scaleDecode(2),5),'yticklabel',num2str(linspace(scaleDecode(1),scaleDecode(2),5)',2),'FontSize',15,'colormap','Jet');
%%
    figureSave
    

