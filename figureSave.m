% saveFigName


if ~exist('saveFigName','var')
    saveFigName = savefilename;
end

if ~exist('Pl','var')
    Pl = PLOT;
end

saveFigName = removePeriodFromName(saveFigName);

if ~exist('D','var')
    D = DIR;
end

if Pl.printPNG
        fileExt = 'png';
        print(gcf,['-d' fileExt],[D.fig saveFigName '.' fileExt])
end

if Pl.printEPS
        fileExt = 'eps';
        print(gcf,['-d' fileExt 'c'],[D.fig saveFigName '.' fileExt])
end

if ~isfield('printFIG',Pl)
    Pl.printFIG=0;
end
if Pl.printFIG
%         fileExt = 'fig';
%         print(gcf,['-d' fileExt],[D.fig saveFigName '.' fileExt])
        fileExt = 'fig';
        savefig([D.fig saveFigName '.' fileExt])
end
