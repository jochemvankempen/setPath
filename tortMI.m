function MI = tortMI(nBins, lowFreqPhase, highFreqPow, plotDistribution)

if nargin < 4
    plotDistribution = 0;
end

%define the edges of the bines
edges = linspace(0,360,nBins+1);

%bin the phases of the lower frequency
[N lowFreqPhaseBin] = histc(mod(lowFreqPhase/pi*180,360),edges);

%compute the composite time series (canolty et al. do this in complex
%plane! difference?) Tort et al. 2010 is unclear in the way they calculate this,
%however in Tort et al. 2008 they just take the composite and comment on
%the way canolty 2006 do it
composite_time_series = highFreqPow.*exp(1i*lowFreqPhase); 

% Create a uniform distribution (U) as a reference to P (below)
U = 1/nBins;

%calculate the mean gamma power in each of the phase bins. Here I take the
%absolute value. Tort et al are unclear about this. If the previous step is
%in complex plane, definitely take absolute value  
% meanGammaPowerBin = zeros(1,length(unique(lowFreqPhaseBin)));

meanGammaPowerBin = zeros(1,nBins);
for ibin = 1:nBins
    if find(ismember(unique(lowFreqPhaseBin),ibin))
        meanGammaPowerBin(ibin) = abs(mean(squeeze(composite_time_series(lowFreqPhaseBin==ibin))));
    else
        meanGammaPowerBin(ibin) = U;
    end
end

% if length(unique(meanGammaPowerBin)) < nBins
%     MI = nan;
%     return
% %     meanGammaPowerBin = zeros(1,nBins);
% end

% Compute Normalized mean amplitude (P(j) in Tort), discrete distribution
P = zeros(1,length(unique(nBins)));
for iBin = 1:nBins
    P(iBin) = meanGammaPowerBin(iBin)/sum(meanGammaPowerBin);
end

%distribution of power over phase bins
if plotDistribution == 1
    figure(1),clf
    bar([P P],'grouped')%, nBins)
    xlim([0 2*length(P)+1])
    set(gca,'XtickLabel',round(linspace(0,2*edges(end),length(get(gca,'xtick'))) ))
end

% Compute the KL (Kullback-Leibler) distance of discrete distribution P
% from a distribution U
D_KL = zeros(1,length(unique(nBins)));
for iBin = 1:length(unique(lowFreqPhaseBin))
    D_KL(iBin) = P(iBin) * log10(P(iBin)/U);
end
D_KL = sum(D_KL);

% Compute the MI
MI = D_KL/log10(nBins);

% calculating D_KL based on shannon distribution of entropy

% for iBin = 1:length(unique(lowFreqPhaseBin))
% H(iBin) = P(iBin) * log10(P(iBin));
% end
% H = - sum(H);
% 
% D_KL = log10(nBins) - H;
% 
% 
% % Compute the MI
% MI = D_KL/log10(nBins);