function [lowFreqFiltThrough highFreqPowThrough, ff] = powerOverPhaseThrough(lowFreqPhase, highFreqPow, lowfreq, fs, windowAroundThrough, plotFigures)
%Help
%This function filters raw data in a lower frequency band for phase values
%and in higher frequency bands (normalized) for power values. It further
%alligns the lower and higher frequency data to the through of the lower
%frequency phase values
%
%INPUT
%data = raw data (1*n)
%lowfreq = variable, filter wil take .5 window around this frequency
%          vector, bandpass between these 2 numbers
%highFreqRange = vector (1*2) indicating the range of the higher
%               frequencies
%highFreqStep = the stepsize of the bandpass filter
%highFreqOverlap = the overlap of the bandpassed higher frequencies
%fs = sampling rate
%windowAroundThrough = window around the through of the lower frequency
%phase in ms
%plotFigures = plot the figures with filtered lower frequency, default = 0
%
%OUTPUT
%lowFreqFiltThrough = Lower Frequencies (amount of througs*(windoAroundThrough*2))
%alligned at the through of this frequency
%highFreqPowThrough = higher frequencies (frequencies *
%(windoAroundThrough*2)), already meaned over amount of throughs
%highFreq2test = 2*n matrix of filtered high frequency bands


if nargin < 5
plotFigures = 0;
end

%% filter for the lower frequency
% if length(lowfreq) == 2
%     lowfreqfilt = lowfreq;
% elseif length(lowfreq) == 1
%     lowfreqfilt = [lowfreq-0.5 lowfreq+0.5];
% else error('lowfreq not specified')
% end

[lowFreqFilt lowFreqPhase] = filtForMI(lowfreqfilt(1),lowfreqfilt(2),data,fs,'phase');

%find the through of the lower frequency, this is -pi
lowFreqThrough = find(lowFreqPhase < (-pi+0.05));

%get a time window around all the lower frequency throughs
lowFreqFiltThrough = zeros(length(lowFreqThrough),2*windowAroundThrough+1);
for ithrough = 1:length(find(lowFreqThrough))
    if lowFreqThrough(ithrough) > windowAroundThrough && lowFreqThrough(end)-lowFreqThrough(ithrough) > windowAroundThrough
        lowFreqThroughTime = lowFreqThrough(ithrough)-windowAroundThrough:1:lowFreqThrough(ithrough)+windowAroundThrough;
        lowFreqFiltThrough(ithrough,:) = lowFreqFilt(lowFreqThroughTime);
    end
end
badThroughs = lowFreqFiltThrough(:,1) == 0;
lowFreqFiltThrough(badThroughs,:) = [];

% gaussian = gauss(windowAroundThrough*2+1,2);

if plotFigures == 1
    figure
%     subplot(3,1,1)
    plot(mean(lowFreqFiltThrough,1))
    xlim([0 windowAroundThrough*2+1])
%     subplot(3,1,2)
%     plot(gaussian(1:(windowAroundThrough*2)))
%     xlim([0 windowAroundThrough*2+1])
%     subplot(3,1,3)
%     plot(mean(lowFreqFiltThrough,1).*gaussian,'r');
%     xlim([0 windowAroundThrough*2+1])
end

%% Filter for the higher frequencies

% if length(highFreqRange) < 2
%     error('highfreq not specified, specify a range (in vector)')
% end
% 
% highFreq = highFreqRange(1):highFreqOverlap:highFreqRange(2)-highFreqStep;
% 
% highFreq2test   = zeros(length(highFreq),2);
% for ihighFreq = 1:length(highFreq)
%     highFreq2test(ihighFreq,:)   = [highFreq(ihighFreq) highFreq(ihighFreq)+highFreqStep];
% end
% 
% highFreqFilt    = zeros(size(highFreq2test,1),length(data));
% highFreqPow     = zeros(size(highFreq2test,1),length(data));
% 
% % parfor ihighFreq = 1:size(highFreq2test,1)
% %     [highFreqFilt(ihighFreq,:) highFreqPow(ihighFreq,:)] = ...
% %         filtForMI(highFreq2test(ihighFreq,1),highFreq2test(ihighFreq,2),data,fs,'power',1);
% % end


movingwin = [0.1 0.001];
params.tapers = [3 5];
params.Fs = 1000;
disp('calculating higher frequencies');
[S,tt,ff] = mtspecgramc(data,movingwin, params);
logS = 10*log10(S);




%%
highFreqPowThrough = zeros(size(ff,2),length(lowFreqThrough),2*windowAroundThrough+1);
for ithrough = 1:length(find(lowFreqThrough))
    if lowFreqThrough(ithrough) > windowAroundThrough && lowFreqThrough(end)-lowFreqThrough(ithrough) > windowAroundThrough
        [~, highFreqThrough] = min(abs(tt*1000-lowFreqThrough(ithrough)));
        
        highFreqThroughTime = highFreqThrough-windowAroundThrough:1:highFreqThrough+windowAroundThrough;
                
        for ihighFreq = 1:size(ff,2)
            highFreqPowThrough(ihighFreq,ithrough,:) = logS(highFreqThroughTime, ihighFreq);
        end
    end
end
badThroughs = highFreqPowThrough(1,:,1) == 0;
highFreqPowThrough(:,badThroughs,:) = [];  

highFreqPowThrough = squeeze(mean(highFreqPowThrough,2));
for ihighFreq = 1:size(ff,2)
    highFreqPowThrough(ihighFreq,:) = squeeze(highFreqPowThrough(ihighFreq,:));%.*gaussian;
end


end



















