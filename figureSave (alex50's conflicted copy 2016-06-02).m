DIR.figPNG = [DIR.fig 'PNG' filesep];
if ~exist([DIR.figPNG],'dir')
    mkdir(DIR.figPNG)
end

DIR.figEPSC = [DIR.fig 'EPS' filesep];
if ~exist([DIR.figEPSC],'dir')
    mkdir(DIR.figEPSC)
end

savefilename = removePeriodFromName(savefilename);

if PLOT.printPNG
        fileExt = 'png';
        print(gcf,['-d' fileExt],[DIR.figPNG savefilename '.' fileExt])
end

if PLOT.printEPS
        fileExt = 'eps';
        print(gcf,['-d' fileExt 'c'],[DIR.figEPSC savefilename '.' fileExt])
end

