function h = plotBrokenVector(xData, yData, idx, plSettings)
%function to plot parts of a vector, e.g. for significance values in a time
%plot.
%
% INPUT
% xData: data on x-axis, e.g. time values
% yData: data on y-axis
% idx: vector of zeros and ones of same length as x and y. Where idx is 1, it will plot yData 

if sum(idx)==0 % if there are no ones, do not plot anything
    h=[];
    return
end

if nargin<4
    plSettings.lineWidth = 2;
    plSettings.color = [0 0 0];
end

idx = idx(:);

if length(yData)==1
    yData = repmat(yData,1,length(xData));
end

changeIdx = find(diff(idx))+1;

if idx(1)==1
    changeIdx = [1; changeIdx];
end
if idx(end)==1
    changeIdx = [changeIdx; length(idx)];
end


for iChange = 1:2:length(changeIdx)-1
    h = plot(xData(changeIdx(iChange):changeIdx(iChange+1)),yData(changeIdx(iChange):changeIdx(iChange+1)),'lineWidth',plSettings.lineWidth,'color',plSettings.color);
end


