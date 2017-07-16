function lp = linearProjection(erp, RT, tt, plot_lp, figureNumber)
% lp = linearProjection(erp, RT, tt, plot_lp)
%
% computes the linear projection of single trials on the average vector.
% 
% INPUT
% erp: time x trial matrix
% RT: reaction time values in samples corresponding to tt (optional). If given, single trial projections
%   will be corrected for their own RT
% tt: time values (optional). If given, performs check of erp matrix
%   dimensions and is used for plotting
% plot_lp: if 1 (optional, default is 0), it will plot the average vector, the norm
%   of the vector, every single trial and the single trial linear projection.
%
% OUTPUT
% lp: linear projection for each individual trial
%
% based on papers 10.1038/78856, 10.1073/pnas.1317557111 and 10.1111/ejn.12859
% jochem van kempen 22/02/2017

if nargin<5 || isempty(figureNumber)
    figureNumber=1;
end
if nargin<4 || isempty(plot_lp)
    plot_lp=0;
end
if nargin<3 || isempty(tt)
    tt= 1:size(erp,1);
else
    if tt(1) ~= 0
        if ~isempty(RT) || nargin<2
            warning('timepoint 1 is not 0, so RTs will not be properly assigned')
            keyboard
        end
    end
end
if nargin<2
    RT=[];
end

% check if erp matrix has correct number of dimension
if ndims(erp)>2
    error('linear projection only works on 1 channel')
end

% check if erp matrix has size (time * trial), if not transpose
if size(erp,1) ~=length(tt)
    disp('transposing erp matrix to fit time x trial')
    erp = erp';
end
[nTime, nTrial] = size(erp);


averageVector   = mean(erp,2); %column vector, average across trials
normVector      = averageVector/sqrt(averageVector'*averageVector)^2;% norm of the average vector
% hold on! post-response data goes into the calculation of the average vector, including for short RT trials (possible contamination of the average vector with artifacts).


% tIdx = (tt>=set.BL(1) & tt<=1500); % this is arbitrary so far.

if plot_lp
    figure(figureNumber),clf
    set(gcf,'Color','White');
    set(gcf,'position',[119    62   989   634])
    set(gcf,'PaperPositionMode','auto')
    subplot(2,1,1)
    boundedline(tt, mean(erp,2), std(erp,[],2))
    set(gca,'fontsize',12)
    title('average vector','fontsize',14)
    xlabel('time (ms)')
    xlim([0 2000])
    
    subplot(2,1,2)
    plot(tt, normVector)
    set(gca,'fontsize',12)
    title('norm of the vector','fontsize',14)
    xlabel('time (ms)')
    YLIM = get(gca,'ylim')
    set(gca,'ylim',YLIM + [-0.1*YLIM(2) 0.1*YLIM(2)])
    xlim([0 2000])
    
    clrs = distinguishable_colors(nTrial);
    
%     saveFigName = [ ];
%     print(gcf,['-d' figureFileType],[paths.pop 'fig' filesep 'p_level' filesep saveFigName '.' figureFileType])


    
    keyboard
%%
end
    
%%
for itrial = 1:nTrial
    
    if ~isempty(RT) 
        tIdx = tt<=RT(itrial);
    else   
        tIdx = tt<=tt(end);
    end

    % for each trial, compute linear projection by multiplying single trial
    % row vector by the norm and deviding by the length
    lp(itrial,1) = (erp(tIdx,itrial)'*normVector(tIdx))/length(normVector(tIdx));

end



%%
if plot_lp
    figure(figureNumber),clf
    set(gcf,'Color','White');
    set(gcf,'position',[119    62   989   634])
    set(gcf,'PaperPositionMode','auto')
    title(['single trial vector & matching linear projection scalar'],'fontsize',14)
    set(gca,'fontsize',12)
    trIdx = 1:20;
    
    %     keyboard
    for itrial = trIdx
        if ~isempty(RT)
            tIdx = tt<=RT(itrial);
        else
            tIdx = tt<=tt(end);
        end
        
        
        %         subplot(2,1,1)
        hold on
        plot(tt(1,tIdx),erp(tIdx,itrial),'color',clrs(itrial,:),'linewidth',1.5)
        %         if ~isempty(RT)
        %             title(['single trial vector, RT=' num2str(RT(itrial))])
        %         else
        %             title(['single trial vector'])
        %         end
        
        %         set(gca,'ylim',YLIM + [-0.1*YLIM(2) 0.1*YLIM(2)])
        xlim([0 1000])
    end
    ylabel('Pupil Diameter','fontsize',14)
    xlabel('Time (ms)','fontsize',14)
    ax1 = gca;

    
    ax2 = axes('Position',get(ax1,'Position'),...
        'XAxisLocation','top',...
        'YAxisLocation','right',...
        'Color','none',...
        'XColor','k','YColor','k','XTickLabel',[],'YTickLabel',[]);
    %         ax2.YLim = [];
    linkaxes([ax1 ax2],'x');
    hold on
%     set(gca,'fontsize',10)
    ylabel('linear projection','fontsize',14)
    for itrial = trIdx
    %         hold on
    plot(900, lp(itrial,1),'.','markersize',20,'color',clrs(itrial,:))
    %         title('single trial linear projection')
    end
%     pause
end

