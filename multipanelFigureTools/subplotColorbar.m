function colH = subplotColorbar(axH,colors,cLimits,location)
% SUBPLOTCOLORBAR Add a legend to a multipanel figure
%   colH = SUBPLOTCOLORBAR(axH,colors,cLimits,location) adds a colorbar in
%   the axes axH.  colors is the desired colormap, with color limits
%   cLimits.  location is 'south', 'east', etc. defined in colorbar.m, and
%   determines the position of the new colorbar within axH.  colH is the
%   handle of the new colorbar.
%
%   Written by Michael Schwendeman, 2014

set(axH,'visible','off')
colormap(colors)
colH = colorbar('peer',axH,'location',location);
set(axH,'cLim',cLimits)