function h = plotBrokenVector(xData, yData, idx, varargin)
% h = plotBrokenVector(xData, yJitter, idx, plSettings)
%
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
    h=NaN;
    return
end

%%% find colors
[found, val, varargin] = parseparam(varargin, 'color');
if found
    colors = val;
else
    colors = [0 0 0];
end

%%% linewidth
[found, val, varargin] = parseparam(varargin, 'linewidth');
if found
    linewidths = val;
else
    linewidths = 1;
end

idx = idx(:);

Y = NaN(length(xData),1);
Y(idx) = yData;

% if exist
h = plot(xData,Y,'lineWidth',linewidths,'color',colors);


% Y=yData;
% if length(Y)==1
% %     yData = repmat(Y,1,length(xData));
%     yData = repmat(Y,1,length(xData));
% end
% 
% Y
% 
% 
% changeIdx = find(diff(idx))+1;
% 
% if idx(1)==1
%     changeIdx = [1; changeIdx];
% end
% if idx(end)==1
%     changeIdx = [changeIdx; length(idx)];
% end
% 
% 
% for iChange = 1:2:length(changeIdx)-1
%     h = plot(xData(changeIdx(iChange):changeIdx(iChange+1)),yData(changeIdx(iChange):changeIdx(iChange+1)),'lineWidth',plSettings.lineWidth,'color',plSettings.color);
% end
% 
% end
% 
% function Y=findMinY(yJitter)
%     % The significance bar needs to be plotted a reasonable distance above all the data points
%     % found over a particular range of X values. So we need to find these data and calculat the 
%     % the minimum y value needed to clear all the plotted data present over this given range of 
%     % x values. 
%     %
%     % This version of the function is a fix from Evan Remington
%     oldXLim = get(gca,'XLim');
%     oldYLim = get(gca,'YLim');
% 
%     axis(gca,'tight')
%     set(gca,'xlim',oldXLim) %Matlab automatically re-tightens y-axis
% 
%     yLim = get(gca,'YLim'); %Now have max y value of all elements within range.
%     Y = min(yLim);
%     
% %     if sign<0
%         Y = Y - abs(0.1 * yJitter * range(yLim));
% %     elseif sign>0
% %         Y = Y + abs(0.1 * range(yLim));
% %     end
%     
%         axis(gca,'normal')
%         set(gca,'XLim',oldXLim,'YLim',[Y - abs(0.2*range(yLim)) oldYLim(2)])
% 
% end %close findMinY



%--------------------
% Parse optional 
% parameters
%--------------------

function [found, val, vars] = parseparam(vars, param)

isvar = cellfun(@(x) ischar(x) && strcmpi(x, param), vars);

if sum(isvar) > 1
    error('Parameters can only be passed once');
end

if any(isvar)
    found = true;
    idx = find(isvar);
    val = vars{idx+1};
    vars([idx idx+1]) = [];
else
    found = false;
    val = [];
end
