function [binEdges, binIdx, binValues] = binVal(inpData, nBin, binType)
% sort values (e.g. RTs) in bins.
%
% INPUT:
% RTs = values to sort
% nBin = number of bins to sort it in
% binType = 'median' or 'equal' split
%
% OUTPUT:
% binEdges = Edges of the bins
% binIdx = idices of bin for every original value
% binValues = original values, binned
%
% Jochem van Kempen, 07/02/2017


[tmpDat,sortIdx] = sort(inpData(:)); % sort all values

[~, revIdx] = sort(sortIdx);%unsort your data

if sum(tmpDat(revIdx) == inpData(:)) ~= length(inpData)
    error('something went wrong sorting')
end      

switch binType
    case 'median'
        if nBin ~= 2
            error('cannot do median split with number of bins unequal to two')
        end
        binEdges = [0 median(tmpDat) tmpDat(end)];
        
    case 'equal'
        binSize    = floor(length(tmpDat)/nBin);
        for iBin = 1:nBin-1
            binEdges(iBin) = tmpDat(binSize*iBin);
        end
        binEdges = [0 binEdges tmpDat(end)];
end

if nargout == 1 % no need for 
    binIdx = [];
    binValues = [];
    return
else
    error('needs some work')
end


binIdx      = zeros(length(inpData),1);
binValues   = cell(nBin,1);

for iBin = 1:nBin
    binIdx(inpData > binEdges(iBin) & inpData <= binEdges(iBin+1))=iBin;
    
    binValues{iBin} = [inpData(binIdx==iBin)];
end


    
    
    
    
    