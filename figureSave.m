% saveFigName


saveFigName = removePeriodFromName(saveFigName);

if Pl.printPNG
        fileExt = 'png';
        print(gcf,['-d' fileExt],[D.fig saveFigName '.' fileExt])
end

if Pl.printEPS
        fileExt = 'eps';
        print(gcf,['-d' fileExt 'c'],[D.fig saveFigName '.' fileExt])
end

if Pl.printFIG
%         fileExt = 'fig';
%         print(gcf,['-d' fileExt],[D.fig saveFigName '.' fileExt])
        fileExt = 'fig';
        savefig([D.fig saveFigName '.' fileExt])
end
