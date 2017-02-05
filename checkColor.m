function checkColor(rgbValues)


yPlot = [0.1*size(rgbValues,1):-0.1:0];

figure(1),clf
hold on
for iCol = 1:size(rgbValues,1)
    plot([0 1],[yPlot(iCol) yPlot(iCol)],'linewidth',20,'color',rgbValues(iCol,:))
    
    
end
ylim([min(yPlot) max(yPlot)])


% Pl.colorblind = [ 
%          0         0    1.0000
%     1.0000         0         0
%     1.0000    1.0000         0
%     0.6602    0.6602    0.6602
%          0         0         0
%     1.0000    0.6445         0
%     1.0000         0    1.0000
%          0    0.5000    0.5000
%          0         0    0.5430
%          0    0.3906         0
%          0    1.0000    1.0000
%     0.5977    0.1953    0.7969];
% 
% Pl.col_ind  = distinguishable_colors(40);
