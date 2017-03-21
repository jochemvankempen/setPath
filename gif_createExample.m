% gifTest

tspan = linspace(0,2*pi,30); tspan(end) = [];
x = linspace(-pi,pi);
figure, j = 0;
M = struct('cdata',{},'colormap',{});
for t = tspan
    j = j+1;
    y1 = sin(x+t); y2 = cos(x-t);
    plot(x,y1,'b',x,y2,'r');
    M(j) = getframe;
end
%
% use cdata from first two frames only
movie2gif(M,{M(1:2).cdata},'test.gif','delaytime',0.05,'loopcount',inf);




t = linspace(0,2*pi,31);
A = sin(t);

iM = 0;
for iT = 1:length(t)
    iM=iM+1;
    
    plot(t(1:iT),A(1:iT))
    xlim([0 2*pi])
    ylim([-1 1])
    M(iM) = getframe;
    
end
    
    
movie2gif(M,{M(1:2).cdata},'sinusTest.gif','delaytime',0.05,'loopcount',inf);
