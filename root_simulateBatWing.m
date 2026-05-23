clc;clear;close all;

% load in a bat wing
load("BatData_Riskin_Cynopterus_02_ds3_sRng6en3_8en3_SpherePurts_angStdv1Deg_setUpForSim.mat");

%% set some params
numMeshPts=length(meshFracMasses);
totalPatMass=.01; % summed mass of all patagia
masses=meshFracMasses*totalPatMass;
% set resting spring lengths
ssfSpan=.6; ssfCord=.8;% spring pretension factors (1+prestrain)
ssfMat=ssfSpan*springDirSpanMat.^2 + ssfCord*springDirCordMat.^2;
L0mat=ssfMat.*refSpringLengthMat;
L0mat(edgeMeshCmat)=L0mat(edgeMeshCmat)/2;% edges more pretensed
% set spring stiffnesses
wbar=mean(refSpringLengthMat,"all","omitmissing");% mean ref spring length, will use as the width of each spring
memThickness=40e-6;
Espan=2.5e4; Ecord=2.5e5; % youngs moduli for span and chord-wise directions
Kspan=Espan*wbar*memThickness;% stiffness of a unit length spring oriented spanwise
Kcord=Ecord*wbar*memThickness;% stiffness of a unit length spring oriented cordwise
Kmat=(Kspan*springDirSpanMat.^2 + Kcord*springDirCordMat.^2)./L0mat;
Kmat(edgeMeshCmat)=Kmat(edgeMeshCmat)*5;% free edges stiffer
% set spring damping factors
Bspan=.0038; Bcord=.0114;% damping for unit length springs oriented span and cordwise
Bmat=(Bspan*springDirSpanMat.^2 + Bcord*springDirCordMat.^2)./L0mat;
% set aerodynamic things
rho=1.225; % air density
A=1.64;% parslew & crowther 2010
B=1.135;% parslew & crowther 2010
C=1.05;% parslew & crowther 2010
mu=.25;% we're saying quarter cord is center of aerodynamic force for blades
a=1.25;
b=a*(1-mu)/mu;% enforces desired aerodynamic center
distribFunc=@(x) betapdf(x,a,b);% using beta function to allocate forces along cord
numInterpSteps=10;
numDataFrames=length(tList);
tList=linspace(0,tList(end),numDataFrames*numInterpSteps); % simulation outputs more steps than the actual bat data
dispTinds=round(linspace(1,length(tList),numDataFrames));
dispTlist=tList(dispTinds);
outputFun=@(t,y,flag)odeProgressReporter(t,y,flag,dispTlist);% for outputing sim calculation progress to command window
%tList=tList(1:200);
numSims=length(meshPoses0List);
meshPoseListCell=cell(numSims,1);
meshVelListCell=cell(numSims,1);
%%
tic;
parfor simNum=1:numSims
    mrkPoseFun=mrkPoseFunList{simNum};
    mrkVelFun=mrkVelFunList{simNum};
    odeFun=@(t,State)batWingODEfun(t,State,mrkPoseFun,...
        mrkVelFun,boneMrkIndPairs,boneWsCell,allBoneMeshInds,allPatMeshInds,Kmat,Bmat,L0mat,...
        elemMembership,bladeElemInfo,A,B,C,rho,distribFunc,masses);
    %% run simulation
    meshPoses0=meshPoses0List{simNum};
    meshVels0=meshVels0List{simNum};
    StateMat0=[meshPoses0(allPatMeshInds,:);meshVels0(allPatMeshInds,:)];
    State0=StateMat0(:);
    if simNum==1
        options=odeset('OutputFcn',outputFun);% just output progress for one, assume they all take about the same amount of time
        [~,stateList] = ode15s(odeFun,tList,State0,options);
    else
        [~,stateList] = ode15s(odeFun,tList,State0);
    end
    
    %% repackage
    % repackage pat mesh states to be slotted into the overall mesh states
    numPatMeshPts=length(allPatMeshInds);
    numBoneMeshPts=length(allBoneMeshInds);
    StateMatList=reshape(stateList',2*numPatMeshPts,3,[]);
    patPoseList=StateMatList(1:numPatMeshPts,:,:);
    patVelList=StateMatList(numPatMeshPts+1:end,:,:);
    % get bone mesh states
    bonePoseList=zeros(numBoneMeshPts,3,size(patPoseList,3));
    boneVelList=bonePoseList;
    for k=1:size(patPoseList,3)
        [boneMeshPoses,boneMeshVels]=getBoneMeshState(tList(k),mrkPoseFun,mrkVelFun,boneMrkIndPairs,boneWsCell);
        bonePoseList(:,:,k)=boneMeshPoses;
        boneVelList(:,:,k)=boneMeshVels;
    end
    % slot everything into the over all mesh state lists
    meshPoseList=zeros(numMeshPts,3,size(patPoseList,3));
    meshVelList=meshPoseList;

    meshPoseList(allBoneMeshInds,:,:)=bonePoseList;
    meshPoseList(allPatMeshInds,:,:)=patPoseList;

    meshVelList(allBoneMeshInds,:,:)=boneVelList;
    meshVelList(allPatMeshInds,:,:)=patVelList;

    meshPoseListCell{simNum}=meshPoseList;
    meshVelListCell{simNum}=meshVelList;

end
runTime=toc;
%% animate one of the runs
simToView=randi(numSims);
meshPoseList=meshPoseListCell{simToView};
bonePoseList=meshPoseList(allBoneMeshInds,:,:);
patPoseList=meshPoseList(allPatMeshInds,:,:);

numFrames=size(meshPoseList,3);
%get plot limits
coordUBs=max(meshPoseList,[],[1 3]);
coordLBs=min(meshPoseList,[],[1 3]);
coordRngs=coordUBs-coordLBs;
maxLims=coordUBs+coordRngs*.2;
minLims=coordLBs-coordRngs*.2;
figure()
Hb=scatter3(bonePoseList(:,1,1),bonePoseList(:,2,1),bonePoseList(:,3,1),"filled");
hold on;
Hp=scatter3(patPoseList(:,1,1),patPoseList(:,2,1),patPoseList(:,3,1),'.');
axis equal
xlim([minLims(1) maxLims(1)]);
ylim([minLims(2) maxLims(2)]);
zlim([minLims(3) maxLims(3)]);
xlabel('x');ylabel('y');zlabel('z');
hold off;


for frame=1:numFrames

    xb=bonePoseList(:,1,frame); yb=bonePoseList(:,2,frame); zb=bonePoseList(:,3,frame);
    xp=patPoseList(:,1,frame); yp=patPoseList(:,2,frame); zp=patPoseList(:,3,frame);
    Hb.XData=xb; Hb.YData=yb; Hb.ZData=zb;
    Hp.XData=xp; Hp.YData=yp; Hp.ZData=zp;

    drawnow
    pause(.01);
end

%% save

save('Riskin_Cynopterus_02_th4en5_Es25e3_Ec25e4_Bs38en4_Bc114en4_SFs6en1_SFc8en1_Edge2xpt5xstiff_rho1225en3_SpherePurts_angStdv1Deg_simResults.mat')

