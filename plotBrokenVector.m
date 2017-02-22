function h = plotBrokenVector(xData, yJitter, idx, plSettings)
%function to plot parts of a vector, e.g. for significance values in a time
%plot.
%
% INPUT
% xData: data on x-axis, e.g. time values
% yData: data on y-axis
% idx: vector of zeros and ones of same length as x and y. Where idx is 1, it will plot yData 
%
% OUTPUT
% h: handle to line

if sum(idx)==0 % if there are no ones, do not plot anything
    h=[];
    return
end

if nargin<4
    plSettings.lineWidth = 2;
    plSettings.color = [0 0 0];
end

idx = idx(:);

Y=findMinY(yJitter);
if length(Y)==1
%     yData = repmat(Y,1,length(xData));
    yData = repmat(Y,1,length(xData));
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

end

function Y=findMinY(yJitter)
    % The significance bar needs to be plotted a reasonable distance above all the data points
    % found over a particular range of X values. So we need to find these data and calculat the 
    % the minimum y value needed to clear all the plotted data present over this given range of 
    % x values. 
    %
    % This version of the function is a fix from Evan Remington
    oldXLim = get(gca,'XLim');
    oldYLim = get(gca,'YLim');

    axis(gca,'tight')
    set(gca,'xlim',oldXLim) %Matlab automatically re-tightens y-axis

    yLim = get(gca,'YLim'); %Now have max y value of all elements within range.
    Y = min(yLim);
    Y = Y - abs(Y*0.1*yJitter);

    
    axis(gca,'normal')
    set(gca,'XLim',oldXLim,'YLim',[Y- abs(Y*0.1*yJitter) oldYLim(2)])

end %close findMinY