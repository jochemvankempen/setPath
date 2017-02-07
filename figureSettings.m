% figureSettings
% scrsz = get(0,'ScreenSize');
% close all

% figure('Position',[1 scrsz(4) scrsz(3)/2 scrsz(4)]),clf
iptsetpref('ImshowBorder','tight');
set(0,'DefaultFigureMenu','none');
format compact;
close all
try
    figure('name','1000','visible', Pl.visible);
catch
    figure
end


% try
%     pos = [ -1549         217        1500         900];
% catch
%     pos = [ 50         150         1500         900];
% end
% set(gcf, 'Position', pos)
set( gcf, 'menubar', 'figure' )
% set(hF,'units','normalized','outerposition',[0.1 0.1 0.6 0.6])


clf
% hFig = figure(1000);
set(gcf,'Color','White');
set(gca,'FontSize',15)
% set(gcf,'Plosition',[1 scrsz(4) scrsz(3) scrsz(4)])

% set(gcf,'DefaultLineLineWidth',2)
% set(gca,'DefaultFontSize',15)

