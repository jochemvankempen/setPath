
if isdir('/gpfs/M2Home/projects/Monash052/jochem/')
    addpath(genpath('/gpfs/M2Home/projects/Monash052/jochem/scripts/remoteScriptsGIT/setPathFunctions/'));
    addpath(genpath('/gpfs/M2Home/projects/Monash052/jochem/human/bin/'));
    addpath(genpath('/gpfs/M2Home/projects/Monash052/jochem/scripts/remoteScriptsGIT/decodeRem/'));
    dataDir =       '/gpfs/M2Home/projects/Monash052/jochem/human/localizer/data/jochemData';
    if matlabpool('size') == 0, matlabpool open 6, end    
elseif isdir('/export/kani/jochem/')
    addpath(genpath('/export/kani/jochem/scripts/setPathFunctions/'));
    addpath(genpath('/export/kani/jochem/human/localizer/bin/'));
    addpath(genpath('/export/kani/jochem/scripts/decodeRem/'));
    dataDir =       '/export/kani/jochem/human/localizer/data/jochemData';
    if matlabpool('size') == 0, matlabpool open, end    
elseif isdir('/Users/Jochem/Documents/UvA/')
    addpath(genpath('/Users/Jochem/Dropbox/scripts/remoteScriptsGIT/setPathFunctions/'));
    addpath(genpath('/Users/Jochem/Documents/UvA/Master/ResearchProject2/bin/'));
    addpath(genpath('/Users/Jochem/Dropbox/scripts/remoteScriptsGIT/decodeRem/'));
    dataDir =       '/Users/Jochem/Documents/UvA/Master/ResearchProject2/data/jochemData/';
    if matlabpool('size') == 0, matlabpool open, end
end
