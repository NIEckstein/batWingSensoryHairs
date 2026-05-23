function [] = plotARefConfig(refPts,cMat,color)
[indsFrom,indsTo]=find(tril(cMat));
x=refPts(:,1)';y=refPts(:,2)';
xPlot=[x(indsFrom);x(indsTo)];
yPlot=[y(indsFrom);y(indsTo)];
figure()
if nargin==2
    plot(xPlot,yPlot);
else
    plot(xPlot, yPlot, color);
end
axis equal;
end
