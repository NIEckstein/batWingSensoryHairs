clc;clear;close all;

% load in a simulation setup data file
load("BatData_Riskin_Cynopterus_02_ds3_sRng6en3_8en3_SpherePurts_angStdv1Deg_setUpForSim_2.mat");
numTrajs=length(mrkPoseFunList);
numMrks=size(refMrkPts,1);
% set up animation parameters
timeDilationFactor=10;% how many times slower should animation be than the actual data
frameRate=40;% desired frames per second for the animation (frames per actual second, not per second of data time)
% set up tList given our desired animation parameters
numSteps=round(frameRate*timeDilationFactor*range(tList));
tList=linspace(tList(1),tList(end),numSteps);
% get marker trajectories sampled at the times in tList
mrkHists=zeros(numMrks,3,numSteps,numTrajs);
for tStep=1:numSteps
    for trajNum=1:numTrajs
        mrkHists(:,:,tStep,trajNum)=mrkPoseFunList{trajNum}(tList(tStep));
    end
end

%% do the animation

colors=hsv(numTrajs);
indsFrom = boneMrkIndPairs(:,1);
indsTo = boneMrkIndPairs(:,2);
FOVbuffer=.05;
xLims=[min(mrkHists(:,1,:,:),[],"all")-FOVbuffer max(mrkHists(:,1,:,:),[],"all")+FOVbuffer-freeStreamVel*tList(end)];
yLims=[min(mrkHists(:,2,:,:),[],"all")-FOVbuffer max(mrkHists(:,2,:,:),[],"all")+FOVbuffer];
zLims=[min(mrkHists(:,3,:,:),[],"all")-FOVbuffer max(mrkHists(:,3,:,:),[],"all")+FOVbuffer];
figure()
for trajNum=1:numTrajs
    mrksNow=mrkHists(:,:,1,trajNum);
    x=mrksNow(:,1)'; y=mrksNow(:,2)'; z=mrksNow(:,3)';
    xPlot=[x(indsFrom);x(indsTo)];
    yPlot=[y(indsFrom);y(indsTo)];
    zPlot=[z(indsFrom);z(indsTo)];
    H(:,trajNum)=plot3(xPlot,yPlot,zPlot,'Color',colors(trajNum,:));
    hold on;
end
grid on;
axis equal;
xlim(xLims); ylim(yLims); zlim(zLims);
for tStep=1:length(tList)
    for trajNum=1:numTrajs
        mrksNow=mrkHists(:,:,tStep,trajNum);
        x=mrksNow(:,1)'; y=mrksNow(:,2)'; z=mrksNow(:,3)';
        xPlot=[x(indsFrom);x(indsTo)];
        yPlot=[y(indsFrom);y(indsTo)];
        zPlot=[z(indsFrom);z(indsTo)];
        for seg=1:size(H,1)
            H(seg,trajNum).XData=xPlot(:,seg);
            H(seg,trajNum).YData=yPlot(:,seg);
            H(seg,trajNum).ZData=zPlot(:,seg);
        end
    end
    xlim(xLims+freeStreamVel*tList(tStep));
    drawnow;
    FlipBook(tStep)=getframe(gcf);

end

% animationFileName='animationOfSkeletonPurturbations';
% Writer=VideoWriter(animationFileName,'MPEG-4');
% 
% 
% Writer.FrameRate = frameRate;
% 
% open(Writer);
% writeVideo(Writer,FlipBook);
% close(Writer);