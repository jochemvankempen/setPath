function [H,P,CI,STATS,auroc] = plotScatter(data,idx,idxLabels,axislabels,Pl,statistics);
% plots scatterplot 
%
% INPUT
% data
%
%
%
% TO DO, adapt function for multiple dimension

if nargin<5 || isempty(Pl)
    setPlotPar
end

if nargin<4 || isempty(axislabels)
    axislabels = {'x','y'};
elseif isnan(axislabels{1})
    axislabels = {'',''};
end

if nargin <2 || isempty(idx)
    idx = ones(size(data,1),1);
end
if nargin <3 || isempty(idxLabels)
    for iLab = 1:length(unique(idx))
        idxLabels{iLab} ={''};
    end
end
nDim = ndims(data);

if ~isfield(Pl,'limits')
    Pl.limits = [min(data(:)) max(data(:))];
    plotLim2compare = [1 10 20 30 50 75 100];
    if Pl.limits(1) < 0;
        for iPlcomp = 1:length(plotLim2compare)
            if Pl.limits(2)<plotLim2compare(iPlcomp)
                Pl.limits = [-plotLim2compare(iPlcomp) plotLim2compare(iPlcomp)];
                break
            end
        end
    else
        for iPlcomp = 1:length(plotLim2compare)
            if Pl.limits(2)<plotLim2compare(iPlcomp)
                Pl.limits = [0 plotLim2compare(iPlcomp)];
                break
            end
        end
    end
end
colors = {'r','b','k','g'};
idx2use = unique(idx);
if nDim==2
%     figureSettings
    hold on
    for iCol = 1:length(idx2use)
        col2use = idx2use(iCol);
        if length(unique(idx))==1
            plot(data(idx==col2use,1), data(idx==col2use,2), Pl.marker,'MarkerSize',Pl.markerSize,'color','k')
        else
            plot(data(idx==col2use,1), data(idx==col2use,2), Pl.marker,'MarkerSize',Pl.markerSize,'color',colors{iCol},'linewidth',2)
        end
    end
    plot(Pl.limits,Pl.limits,'k','lineWidth',Pl.lineWidth)
%     lsline
    xlim(Pl.limits)
    ylim(Pl.limits)

    set(gca,'fontsize',Pl.fontsizeAxes)
    xlabel(axislabels{1},'interpret','none','fontsize',Pl.fontsizeLabel)
    ylabel(axislabels{2},'interpret','none','fontsize',Pl.fontsizeLabel)
    axis square
    
    clear H P CI STATS
    for iIdx = 1:length(idx2use)
        auroc(iIdx) = f_auroc(data(idx==iIdx,1),data(idx==iIdx,2));
        switch statistics
            case 'paired'
                [H{iIdx},P{iIdx},CI{iIdx},STATS{iIdx}] = ttest(data(idx==iIdx,1),data(idx==iIdx,2));
            case 'unpaired'
                [H{iIdx},P{iIdx},CI{iIdx},STATS{iIdx}] = ttest2(data(idx==iIdx,1),data(idx==iIdx,2));                
        end
        %             text(mean(Pl.limits),Pl.limits(1) + max(abs(Pl.limits))*(.12*iIdx),sprintf('auroc = %1.3f, (p = %0.3f)',  auroc(iIdx) , P{iIdx}),'color',colors{iIdx})
        if 0
            if length(unique(idx))==1
                text(mean(Pl.limits),Pl.limits(1) + max(abs(Pl.limits))*(.12*iIdx),sprintf('p = %0.3f', P{iIdx}),'color','k')
            else
                text(mean(Pl.limits),Pl.limits(1) + max(abs(Pl.limits))*(.12*iIdx),sprintf('p = %0.3f', P{iIdx}),'color',colors{iIdx})
            end
        else
            text(mean(Pl.limits),Pl.limits(1) + iIdx,[idxLabels{iIdx} ' ' sprintf('p = %0.3f', P{iIdx})],'color',colors{iIdx},'fontsize',Pl.fontsizeAxes)
            
        end
    end
else
    disp('function not ready for more dimensions')
end

% if nDim==3