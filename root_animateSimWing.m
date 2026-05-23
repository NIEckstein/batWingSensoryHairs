clc;clear;close all;

load("Riskin_Cynopterus_02_th4en5_Es25e3_Ec25e4_Bs38en4_Bc114en4_SFs6en1_SFc8en1_Edge2xpt5xstiff_rho1225en3_SpherePurts_angStdv1Deg_simResults.mat");
%% animate one of the runs
% simToView=randi(numSims);% pick a run to view at random
 simToView=3;% pick a specific run to view
meshPoseList=meshPoseListCell{simToView};
bonePoseList=meshPoseList(allBoneMeshInds,:,:);
patPoseList=meshPoseList(allPatMeshInds,:,:);

numFrames=size(meshPoseList,3);
%get starting plot limits (the viewing window will move with the bat)
coordUBs=max(meshPoseList,[],[1 3]);
% subtract the max free stream translation from the starting x coord UB so our window
% isnt problematically elongated
coordUBs(1)=coordUBs(1)-freeStreamVel*tList(end);
coordLBs=min(meshPoseList,[],[1 3]);
coordRngs=coordUBs-coordLBs;
maxLims0=coordUBs+coordRngs*.2;
minLims0=coordLBs-coordRngs*.2;
figure()
Hb=scatter3(bonePoseList(:,1,1),bonePoseList(:,2,1),bonePoseList(:,3,1),'k',"filled");
hold on;
numBlades=length(bladeElemInfo);
colors=hsv(numBlades);
for bladeNum=1:numBlades
    inds=bladeElemInfo(bladeNum).elemInds;
    mem=elemMembership(inds,:);
    Hp(bladeNum)=patch('Faces',mem,'Vertices',meshPoseList(:,:,1),'FaceColor',colors(bladeNum,:));
end
axis equal
xlim([minLims0(1) maxLims0(1)]);
ylim([minLims0(2) maxLims0(2)]);
zlim([minLims0(3) maxLims0(3)]);
xlabel('x');ylabel('y');zlabel('z');
view(-25,60)
hold off;
% info for the saved mp4 
timeDilationFactor=10;% how many times slower should animation be than the actual data
frameRate=40;% desired frames per second for the animation (frames per actual second, not per second of data time)
desiredNumFramesInMP4=range(tList)*frameRate*timeDilationFactor;
skipNum=round(numFrames/desiredNumFramesInMP4);
dispFrame=0;
freeStreamFrameTrans=-freeStreamVel*tList;
for frame=1:skipNum:numFrames
    % translate viewing window
    freeStreamX=freeStreamFrameTrans(frame);
    xlim([minLims0(1) maxLims0(1)]-freeStreamX);
    % update plot data
    xb=bonePoseList(:,1,frame); 
    yb=bonePoseList(:,2,frame); zb=bonePoseList(:,3,frame);
    Hb.XData=xb; Hb.YData=yb; Hb.ZData=zb;
    Vp=meshPoseList(:,:,frame);
    for bladeNum=1:numBlades
    Hp(bladeNum).Vertices=Vp;
    end
    dispFrame=dispFrame+1;
    drawnow
    FlipBook(dispFrame)=getframe(gcf);
    pause(.05);
end

% animationFileName='animationOfSimulatedWingWithBlades';
% Writer=VideoWriter(animationFileName,'MPEG-4');
% 
% 
% Writer.FrameRate = frameRate;
% 
% open(Writer);
% writeVideo(Writer,FlipBook);
% close(Writer);