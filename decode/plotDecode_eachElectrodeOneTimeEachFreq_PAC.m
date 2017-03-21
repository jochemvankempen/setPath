function  plotDecode_eachElectrodeOneTimeEachFreq_PAC(loadfilename,savefilename,DP,PLOT,DIR,trialTypes)
if ~exist([DIR.decodeCheck filesep loadfilename '.mat'],'file'), return, end

clear hl

load([DIR.decodeCheck filesep loadfilename '.mat'])

if isstruct(times2downsample)
    times2downsample = times2downsample.phase;
end
if isempty(dir(DIR.fig))
    mkdir(DIR.fig)
end

%%
figureSettings
colors = colormap(hsv(4));

if DP.permute || DP.permutePAC
    vDim = [1 2 3 4];
else
    vDim = [1 2 3];
end


matrix2plot = [pMeanTestCorr(2).TF; NaN(1,length(pMeanTestCorr(2).TF)); pMeanTestCorr(3).TF];
matrix2plot = [[NaN NaN pMeanTestCorr(1).TF]' NaN(1,length(pMeanTestCorr(1).TF)+2)' matrix2plot];
figureSettings

colormap default
imagesc(matrix2plot')
axis xy
colorbar

XTICK = 3:2:length(freq.phase)+2;
YTICK = 3:2:length(freq.power)+2;

XLABEL = [NaN NaN round(freq.phase)];
YLABEL = [NaN NaN round(freq.power)];
set(gca,'XTick',XTICK,'XTickLabel',XLABEL(XTICK))
set(gca,'YTick',YTICK,'YTickLabel',YLABEL(YTICK))

xlabel('phase freq')
ylabel('power freq')

figureSave
