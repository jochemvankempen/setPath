function plotSignificance_bar(YLIM,TEST,cond2compare,Pl)


if ~TEST.H
    return
end

if nargin<4
    Pl.lineLength   = 0.02;
    Pl.markerSize   = 15;
    Pl.lineWidth    = 2;
end

% YLIM = get(gca,'ylim');
YLIM2use = max(YLIM) + (max(YLIM)*0.05);


hold on
plot([cond2compare(1) cond2compare(2)],[max(YLIM2use) max(YLIM2use)],'k','linewidth',Pl.lineWidth)
% plot([cond2compare(1) cond2compare(1)],[max(YLIM2use)-(max(YLIM2use)*Pl.lineLength) max(YLIM2use)+(max(YLIM2use)*Pl.lineLength)],'k','linewidth',Pl.lineWidth)
% plot([cond2compare(2) cond2compare(2)],[max(YLIM2use)-(max(YLIM2use)*Pl.lineLength) max(YLIM2use)+(max(YLIM2use)*Pl.lineLength)],'k','linewidth',Pl.lineWidth)

if TEST.P < 0.001
    signMark = '***';
    markCorrection = 0.2+ 0.2 * (8/Pl.markerSize);
elseif TEST.P < 0.01    
    signMark = '**';
    markCorrection = 0.15+ 0.15 * (8/Pl.markerSize);
elseif TEST.P < 0.05    
    signMark = '*';
    markCorrection = 0.05+ 0.15 * (8/Pl.markerSize);
end

% text(mean(cond2compare)-markCorrection,max(YLIM2use)+(max(YLIM2use)*Pl.lineLength),signMark,'fontsize',Pl.markerSize)
text(mean(cond2compare)-markCorrection,max(YLIM2use)+(max(YLIM2use)*0.005),signMark,'fontsize',Pl.markerSize)



% keyboard

