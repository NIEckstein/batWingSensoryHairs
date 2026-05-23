function [] = plotARefConfig_noNewFig(refPts,cMat,color,lineWidth)
[indsFrom,indsTo]=find(tril(cMat));
x=refPts(:,1)';y=refPts(:,2)';
xPlot=[x(indsFrom);x(indsTo)];
yPlot=[y(indsFrom);y(indsTo)];

if nargin==2
    plot(xPlot,yPlot);
elseif nargin==3
    plot(xPlot, yPlot, color);
else
    plot(xPlot, yPlot, color,'LineWidth',lineWidth);
end
axis equal;
grid on;
end
