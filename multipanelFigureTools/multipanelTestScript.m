%%%%    Test script for multi-panel figure tools    %%%%
%   Five methods for making multi-panel figures with legends and colorbars
%   Written by Michael Schwendeman, 2014
clc
clearvars
close all

% setup test data
xPlot = 0:0.1:1;
stDev = 0:0.02:0.1;
ex = 1:3;
nx = length(xPlot);
ns = length(stDev);
ne = length(ex);
yPlot = zeros(nx,ne,ns);
for i = 1:ne
    for j = 1:ns
        yPlot(:,i,j) = 1+xPlot.^ex(i)+stDev(j)*randn(1,nx);
    end
end

% setup colors
minColor = min(stDev);
maxColor = max(stDev);
useColormap = cool;
lineColors = zeros(ns,3);
for j=1:ns
    lineColors(j,:) = interp2(1:3,linspace(minColor,maxColor,64)',useColormap,1:3,stDev(j),'linear'); 
end

% setup linestyles
lineStyles = {'-','--','-.'};

% setup legend string
legStr = {'Linear','Quadratic','Cubic'};
%% Method 1 - subplot with individual legends and colorbars, Matlab defaults
% individual subplots are crowded and margins are unnecessarily large

figure(1)
figSize = [1 1 8 6];
set(gcf,'units','inches','position',figSize,'paperposition',figSize);

subplot(2,2,1)
plotH = zeros(ne,ns);
hold on
for i = 1:ne
    for j = 1:ns
        plotH(i,j) = plot(xPlot,yPlot(:,i,j));
        set(plotH(i,j),'color',lineColors(j,:),'linestyle',char(lineStyles(i)))
    end
end
hold off
ylabel('y variable')
title('Linear Plot')
legend(plotH(1:ne,1),legStr,'location','southeast')
colormap(useColormap)
colH = colorbar('location','northoutside');
xlabel(colH,'Noise Standard Deviation')
set(gca,'cLim',[minColor maxColor])

subplot(2,2,2)
plotH = zeros(ne,ns);
hold on
for i = 1:ne
    for j = 1:ns
        plotH(i,j) = plot(xPlot,yPlot(:,i,j));
        set(plotH(i,j),'color',lineColors(j,:),'linestyle',char(lineStyles(i)))
    end
end
hold off
title('Lin-Log Plot')
legend(plotH(1:ne,1),legStr,'location','southeast')
colormap(useColormap)
colH = colorbar('location','northoutside');
xlabel(colH,'Noise Standard Deviation')
set(gca,'cLim',[minColor maxColor],'xScale','log')

subplot(2,2,3)
plotH = zeros(ne,ns);
hold on
for i = 1:ne
    for j = 1:ns
        plotH(i,j) = plot(xPlot,yPlot(:,i,j));
        set(plotH(i,j),'color',lineColors(j,:),'linestyle',char(lineStyles(i)))
    end
end
hold off
xlabel('x variable')
ylabel('y variable')
title('Log-Lin Plot')
legend(plotH(1:ne,1),legStr,'location','southeast')
colormap(useColormap)
colH = colorbar('location','northoutside');
xlabel(colH,'Noise Standard Deviation')
set(gca,'cLim',[minColor maxColor],'yScale','log')

subplot(2,2,4)
plotH = zeros(ne,ns);
hold on
for i = 1:ne
    for j = 1:ns
        plotH(i,j) = plot(xPlot,yPlot(:,i,j));
        set(plotH(i,j),'color',lineColors(j,:),'linestyle',char(lineStyles(i)))
    end
end
hold off
xlabel('x variable')
title('Log-Log Plot')
legend(plotH(1:ne,1),legStr,'location','southeast')
colormap(useColormap)
colH = colorbar('location','northoutside');
xlabel(colH,'Noise Standard Deviation')
set(gca,'cLim',[minColor maxColor],'xScale','log','yScale','log')

%% Method 2 - multipanel tools, simplest 
% better, but note space above colorbar, tall legends

figure(2)
figSize = [1 1 8 6];
set(gcf,'units','inches','position',figSize,'paperposition',figSize);
numAxRows = 3;
numAxCols = 3;
gridBorder = [0.6 0 0.6 0];
axBorder = [0 0.6 0 0.6];
axH = nan(numAxRows,numAxCols);

axH(2,1) = smartsubplot(gcf,numAxRows,numAxCols,2,1,gridBorder,axBorder);
plotH = zeros(ne,ns);
hold on
for i = 1:ne
    for j = 1:ns
        plotH(i,j) = plot(xPlot,yPlot(:,i,j));
        set(plotH(i,j),'color',lineColors(j,:),'linestyle',char(lineStyles(i)))
    end
end
hold off
ylabel('y variable')
title('Linear Plot')
set(gca,'xScale','lin','yScale','lin')

axH(2,2) = smartsubplot(gcf,numAxRows,numAxCols,2,2,gridBorder,axBorder);
plotH = zeros(ne,ns);
hold on
for i = 1:ne
    for j = 1:ns
        plotH(i,j) = plot(xPlot,yPlot(:,i,j));
        set(plotH(i,j),'color',lineColors(j,:),'linestyle',char(lineStyles(i)))
    end
end
hold off
title('Lin-Log Plot')
set(gca,'xScale','log','yScale','lin')

axH(1,1) = smartsubplot(gcf,numAxRows,numAxCols,1,1,gridBorder,axBorder);
plotH = zeros(ne,ns);
hold on
for i = 1:ne
    for j = 1:ns
        plotH(i,j) = plot(xPlot,yPlot(:,i,j));
        set(plotH(i,j),'color',lineColors(j,:),'linestyle',char(lineStyles(i)))
    end
end
hold off
xlabel('x variable')
ylabel('y variable')
title('Log-Lin Plot')
set(gca,'xScale','lin','yScale','log')

axH(1,2) = smartsubplot(gcf,numAxRows,numAxCols,1,2,gridBorder,axBorder);
plotH = zeros(ne,ns);
hold on
for i = 1:ne
    for j = 1:ns
        plotH(i,j) = plot(xPlot,yPlot(:,i,j));
        set(plotH(i,j),'color',lineColors(j,:),'linestyle',char(lineStyles(i)))
    end
end
hold off
xlabel('x variable')
title('Log-Log Plot')
set(gca,'xScale','log','yScale','log')

legH = zeros(2,1);
axH(1,3) = smartsubplot(gcf,numAxRows,numAxCols,1,3,gridBorder,axBorder);
legH(1) = subplotLegend(axH(1,3),plotH(1:ne),legStr,'vertical');
axH(2,3) = smartsubplot(gcf,numAxRows,numAxCols,2,3,gridBorder,axBorder);
legH(2) = subplotLegend(axH(2,3),plotH(1:ne),legStr,'vertical');

axH(3,1) = smartsubplot(gcf,numAxRows,numAxCols,3,1:2,gridBorder,axBorder);
colH = subplotColorbar(axH(3,1),useColormap,[minColor maxColor],'south');
xlabel(colH,'Noise Standard Deviation')

%% Method 3 - multipanel tools, slightly more involved 
% Fix gaps by adjusting gridBorder and using separate axBorder for legends

figure(3)
figSize = [1 1 8 6];
set(gcf,'units','inches','position',figSize,'paperposition',figSize);
numAxRows = 3;
numAxCols = 3;
gridBorder = [0.6 -0.6 0.6 -1.2];
axBorder = [0 0.6 0 0.6];
axBorderLeg = [-0.3 0.9 0.6 0.9];
axH = nan(numAxRows,numAxCols);

axH(2,1) = smartsubplot(gcf,numAxRows,numAxCols,2,1,gridBorder,axBorder);
plotH = zeros(ne,ns);
hold on
for i = 1:ne
    for j = 1:ns
        plotH(i,j) = plot(xPlot,yPlot(:,i,j));
        set(plotH(i,j),'color',lineColors(j,:),'linestyle',char(lineStyles(i)))
    end
end
hold off
ylabel('y variable')
title('Linear Plot')
set(gca,'xScale','lin','yScale','lin')

axH(2,2) = smartsubplot(gcf,numAxRows,numAxCols,2,2,gridBorder,axBorder);
plotH = zeros(ne,ns);
hold on
for i = 1:ne
    for j = 1:ns
        plotH(i,j) = plot(xPlot,yPlot(:,i,j));
        set(plotH(i,j),'color',lineColors(j,:),'linestyle',char(lineStyles(i)))
    end
end
hold off
title('Lin-Log Plot')
set(gca,'xScale','log','yScale','lin')

axH(1,1) = smartsubplot(gcf,numAxRows,numAxCols,1,1,gridBorder,axBorder);
plotH = zeros(ne,ns);
hold on
for i = 1:ne
    for j = 1:ns
        plotH(i,j) = plot(xPlot,yPlot(:,i,j));
        set(plotH(i,j),'color',lineColors(j,:),'linestyle',char(lineStyles(i)))
    end
end
hold off
xlabel('x variable')
ylabel('y variable')
title('Log-Lin Plot')
set(gca,'xScale','lin','yScale','log')

axH(1,2) = smartsubplot(gcf,numAxRows,numAxCols,1,2,gridBorder,axBorder);
plotH = zeros(ne,ns);
hold on
for i = 1:ne
    for j = 1:ns
        plotH(i,j) = plot(xPlot,yPlot(:,i,j));
        set(plotH(i,j),'color',lineColors(j,:),'linestyle',char(lineStyles(i)))
    end
end
hold off
xlabel('x variable')
title('Log-Log Plot')
set(gca,'xScale','log','yScale','log')

legH = zeros(2,1);
axH(1,3) = smartsubplot(gcf,numAxRows,numAxCols,1,3,gridBorder,axBorderLeg);
legH(1) = subplotLegend(axH(1,3),plotH(1:ne),legStr,'vertical');
axH(2,3) = smartsubplot(gcf,numAxRows,numAxCols,2,3,gridBorder,axBorderLeg);
legH(2) = subplotLegend(axH(2,3),plotH(1:ne),legStr,'vertical');

axH(3,1) = smartsubplot(gcf,numAxRows,numAxCols,3,1:2,gridBorder,axBorder);
colH = subplotColorbar(axH(3,1),useColormap,[minColor maxColor],'south');
xlabel(colH,'Noise Standard Deviation')

%% Method 4 - multipanel tools, slightly more involved 
% Fix gaps by adjusting numAxRows and numAxCols, have plots extend over
% multiple rows and columns, while legend and colorbar are only one

figure(4)
figSize = [1 1 8 6];
set(gcf,'units','inches','position',figSize,'paperposition',figSize);
numAxRows = 7;
numAxCols = 5;
gridBorder = [0.6 0 0.6 0];
axBorder = [0 0.6 0 0.6];
axBorderLeg = [-0.3 0.3 0.3 0.3];
axH = nan(numAxRows,numAxCols);

axH(2,1) = smartsubplot(gcf,numAxRows,numAxCols,4:6,1:2,gridBorder,axBorder);
plotH = zeros(ne,ns);
hold on
for i = 1:ne
    for j = 1:ns
        plotH(i,j) = plot(xPlot,yPlot(:,i,j));
        set(plotH(i,j),'color',lineColors(j,:),'linestyle',char(lineStyles(i)))
    end
end
hold off
ylabel('y variable')
title('Linear Plot')
set(gca,'xScale','lin','yScale','lin')

axH(2,2) = smartsubplot(gcf,numAxRows,numAxCols,4:6,3:4,gridBorder,axBorder);
plotH = zeros(ne,ns);
hold on
for i = 1:ne
    for j = 1:ns
        plotH(i,j) = plot(xPlot,yPlot(:,i,j));
        set(plotH(i,j),'color',lineColors(j,:),'linestyle',char(lineStyles(i)))
    end
end
hold off
title('Lin-Log Plot')
set(gca,'xScale','log','yScale','lin')

axH(1,1) = smartsubplot(gcf,numAxRows,numAxCols,1:3,1:2,gridBorder,axBorder);
plotH = zeros(ne,ns);
hold on
for i = 1:ne
    for j = 1:ns
        plotH(i,j) = plot(xPlot,yPlot(:,i,j));
        set(plotH(i,j),'color',lineColors(j,:),'linestyle',char(lineStyles(i)))
    end
end
hold off
xlabel('x variable')
ylabel('y variable')
title('Log-Lin Plot')
set(gca,'xScale','lin','yScale','log')

axH(1,2) = smartsubplot(gcf,numAxRows,numAxCols,1:3,3:4,gridBorder,axBorder);
plotH = zeros(ne,ns);
hold on
for i = 1:ne
    for j = 1:ns
        plotH(i,j) = plot(xPlot,yPlot(:,i,j));
        set(plotH(i,j),'color',lineColors(j,:),'linestyle',char(lineStyles(i)))
    end
end
hold off
xlabel('x variable')
title('Log-Log Plot')
set(gca,'xScale','log','yScale','log')

legH = zeros(1,1);
axH(1,3) = smartsubplot(gcf,numAxRows,numAxCols,3:4,5,gridBorder,axBorderLeg);
legH(1) = subplotLegend(axH(1,3),plotH(1:ne),legStr,'vertical');

axH(3,1) = smartsubplot(gcf,numAxRows,numAxCols,7,1:4,gridBorder,axBorder);
colH = subplotColorbar(axH(3,1),useColormap,[minColor maxColor],'south');
xlabel(colH,'Noise Standard Deviation')

%% Method 5 - multipanel tools, most involved, best technique
% individual gridBorder and axBorder for plots, legend, and colorbar
% sections

figure(5)
figSize = [1 1 8 6];
set(gcf,'units','inches','position',figSize,'paperposition',figSize);
numAxRows = 2;
numAxCols = 2;
gridBorderPlots = [0 2 0 1.2];
axBorderPlots = [0.6 0 0.6 0];
gridBorderLeg = [figSize(3)-gridBorderPlots(2) 0 gridBorderPlots(3) gridBorderPlots(4)];
axBorderLeg = [0.3 0.3 axBorderPlots(3)+1.6 axBorderPlots(4)+1.6];
gridBorderCol = [gridBorderPlots(1) gridBorderPlots(2) figSize(4)-gridBorderPlots(4) 0];
axBorderCol = [axBorderPlots(1) axBorderPlots(2) 0.4 0.2];


axH = nan(numAxRows,numAxCols);

axH(2,1) = smartsubplot(gcf,numAxRows,numAxCols,2,1,gridBorderPlots,axBorderPlots);
plotH = zeros(ne,ns);
hold on
for i = 1:ne
    for j = 1:ns
        plotH(i,j) = plot(xPlot,yPlot(:,i,j));
        set(plotH(i,j),'color',lineColors(j,:),'linestyle',char(lineStyles(i)))
    end
end
hold off
ylabel('y variable')
title('Linear Plot')
set(gca,'xScale','lin','yScale','lin')

axH(2,2) = smartsubplot(gcf,numAxRows,numAxCols,2,2,gridBorderPlots,axBorderPlots);
plotH = zeros(ne,ns);
hold on
for i = 1:ne
    for j = 1:ns
        plotH(i,j) = plot(xPlot,yPlot(:,i,j));
        set(plotH(i,j),'color',lineColors(j,:),'linestyle',char(lineStyles(i)))
    end
end
hold off
title('Lin-Log Plot')
set(gca,'xScale','log','yScale','lin')

axH(1,1) = smartsubplot(gcf,numAxRows,numAxCols,1,1,gridBorderPlots,axBorderPlots);
plotH = zeros(ne,ns);
hold on
for i = 1:ne
    for j = 1:ns
        plotH(i,j) = plot(xPlot,yPlot(:,i,j));
        set(plotH(i,j),'color',lineColors(j,:),'linestyle',char(lineStyles(i)))
    end
end
hold off
xlabel('x variable')
ylabel('y variable')
title('Log-Lin Plot')
set(gca,'xScale','lin','yScale','log')

axH(1,2) = smartsubplot(gcf,numAxRows,numAxCols,1,2,gridBorderPlots,axBorderPlots);
plotH = zeros(ne,ns);
hold on
for i = 1:ne
    for j = 1:ns
        plotH(i,j) = plot(xPlot,yPlot(:,i,j));
        set(plotH(i,j),'color',lineColors(j,:),'linestyle',char(lineStyles(i)))
    end
end
hold off
xlabel('x variable')
title('Log-Log Plot')
set(gca,'xScale','log','yScale','log')

legH = zeros(1,1);
axH(1,3) = smartsubplot(gcf,1,1,1,1,gridBorderLeg,axBorderLeg);
legH(1) = subplotLegend(axH(1,3),plotH(1:ne),legStr,'vertical');

axH(3,1) = smartsubplot(gcf,1,1,1,1,gridBorderCol,axBorderCol);
colH = subplotColorbar(axH(3,1),useColormap,[minColor maxColor],'south');
xlabel(colH,'Noise Standard Deviation')

print(5,'-dpng','-r75','/Users/mike/Documents/MATLAB/multipanelFigureTools/finalImage.png')