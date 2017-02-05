function axH = smartsubplot(figH,numAxRows,numAxCols,thisAxRow,thisAxCol,gridBorder,axBorder)
% SMARTSUBPLOT Create axes in tiled positions, with explicit borders
%   axH =
%   SMARTSUBPLOT(figH,numAxRows,numAxCols,thisAxRow,thisAxCol,gridBorder,axBorder)
%   draws an axes in figure figH at a panel defined by thisAxRow and
%   thisAxCol, in a grid defined by numAxRows, numAxCols, gridBorder and
%   axBorder.  thisAxRow and thisAxCol may be scalars, indicating an axes
%   in a single grid cell, or vectors, indicating an axes spanning multiple
%   grid cells. thisAxRow and thisAxCol are defined relative to the
%   bottom-left cell. gridBorder and axBorder are four-element vectors
%   defining the margins of the grid within the figure and the axes within
%   the panel, and are of the form [left, right, bottom, top], all in the
%   current units of the figure. Returns the axes handle, axH.  
%   See schematic below:
%   
%                     numAxCols = 2
%          <------------------------------->
%    _______________________________________________
%   |                    ^                         |
%   |                    | gridBorder(4)           |                                  
%   |      ______________v__________________  gridBorder(2)
%   |     |(2,1) ^         |(2,2)           |<---->|        ^
%   |     |      |axBorder(4)               |      |        |
%   |     |    __V_____    |    ________    |      |        |
%   |     |   |        |<->|   |        |   |      |        |
%   |     |   |        |axBorder(2)     |   |      |        |
%   |     |<->|________|   |   |________|   |      |        |   
%   |  axBorder(1)  ^      |                |      |        |
%   |     |         |axBorder(3)            |      |        | numAxRows = 2
%   |     |_________V______|________________|      |        |
%   |     |thisAxRow = 1   |(1,2)           |      |        |
%   |     |thisAxCol = 1   |                |      |        |         
%   |     |    ________    |    ________    |      |        | 
%   |     |   |        |   |   |        |   |      |        | 
%   |     |   |        |   |   |        |   |      |        |   
%   |     |   |________|   |   |________|   |      |        | 
%   |     |                |                |      |        |
%   |<--->|________________|________________|      |        v
%   | gridBorder(1)          ^                     |                   
%   |                        | gridBorder(3)       |   
%   |________________________v_____________________|  
%
%   Written by Michael Schwendeman, 2014

minRow = min(thisAxRow);
maxRow = max(thisAxRow);
minCol = min(thisAxCol);
maxCol = max(thisAxCol);

figPos = get(figH,'position');
figUnits = get(figH,'units');
figWidth = figPos(3);
figHeight = figPos(4);

gridWidth = figWidth-(gridBorder(1)+gridBorder(2));
gridHeight = figHeight-(gridBorder(3)+gridBorder(4));

% check that grid dimensions are not negative
if gridWidth < 0
    error('smartsubplot:gridWidth',['The grid borders are too large for the current figure width.\n',...
        'The sum of left and right grid borders must be less than the figure width of: %.2f %s.\n',...
        'Decrease the grid borders or increase the figure size to fix the problem.\n'],...
        figWidth,figUnits)
end
if gridHeight < 0
    error('smartsubplot:gridHeight',['The grid borders are too large for the current figure height.\n',...
        'The sum of top and bottom grid borders must be less than the figure height of: %.2f %s.\n',...
        'Decrease the grid borders or increase the figure size to fix the problem.'],...
        figHeight,figUnits)
end

cellWidth = gridWidth/numAxCols;
cellHeight = gridHeight/numAxRows;


cellLeftCorner = gridBorder(1)+(minCol-1)*cellWidth;
cellBottomCorner = gridBorder(3)+(minRow-1)*cellHeight;
cellRightCorner = gridBorder(1)+maxCol*cellWidth;
cellTopCorner = gridBorder(3)+maxRow*cellHeight;

axLeftCorner = cellLeftCorner + axBorder(1);
axWidth = cellRightCorner - axBorder(2) - axLeftCorner;
axBottomCorner = cellBottomCorner + axBorder(3);
axHeight = cellTopCorner - axBorder(4) - axBottomCorner;

% check that axis dimensions are not negative
if axWidth < 0
    error('smartsubplot:axWidth',['The axis borders are too large for the current cell width.\n',...
        'The sum of left and right axis borders must be less than the cell width of: %.2f %s.\n',...
        'Decrease the axis borders, grid borders, or the number of cells, or increase the figure size, to fix the problem.'],...
        cellRightCorner-cellLeftCorner,figUnits)
end
if axHeight < 0
    error('smartsubplot:axHeight',['The axis borders are too large for the current cell height.\n',...
        'The sum of top and bottom axis borders must be less than the cell height of: %.2f %s.\n',...
        'Decrease the axis borders, grid borders, or the number of cells, or increase the figure size, to fix the problem.'],...
        cellTopCorner-cellBottomCorner,figUnits)
end

axH = axes('units',figUnits,'position',[axLeftCorner,axBottomCorner,axWidth,axHeight]);

