function [dStatedt] = batWingODEfun(t,State,mrkPoseFun,mrkVelFun,...
    boneMrkIndPairs,boneWsCell,boneMeshInds,patMeshInds,Kmat,Bmat,L0mat,...
    elemMembership,bladeElemInfo,A,B,C,rho,distribFunc,masses)
numPatMeshPts=length(patMeshInds);
numBoneMeshPts=length(boneMeshInds);
% unpack pat mesh state 
StateMat=reshape(State,[2*numPatMeshPts 3]);
patMeshPoses=StateMat(1:numPatMeshPts,:); 
patMeshVels=StateMat(numPatMeshPts+1:end,:);
% get bone mesh state using the real marker kinematic data
[boneMeshPoses,boneMeshVels]=getBoneMeshState(t,mrkPoseFun,mrkVelFun,...
    boneMrkIndPairs,boneWsCell);
% build total mesh state
meshPoses=zeros(numBoneMeshPts+numPatMeshPts,3);
meshVels=meshPoses;

meshPoses(boneMeshInds,:)=boneMeshPoses;
meshVels(boneMeshInds,:)=boneMeshVels;

meshPoses(patMeshInds,:)=patMeshPoses;
meshVels(patMeshInds,:)=patMeshVels;

% compute the forces on the mesh masses
mechForces=getMechForces(meshPoses,meshVels,Kmat,Bmat,L0mat);
airForces=getAirForces(meshPoses,meshVels,elemMembership,bladeElemInfo,A,B,C,rho,distribFunc);
Forces=mechForces+airForces;

% compute the mesh accelerations
meshAccs=Forces./masses;

patMeshAccs=meshAccs(patMeshInds,:);

% repackage
dStateMatdt=[patMeshVels;patMeshAccs];
dStatedt=dStateMatdt(:);

% % add in some gravity;
% meshAccs(:,3)=meshAccs(:,3)-9.81;


end