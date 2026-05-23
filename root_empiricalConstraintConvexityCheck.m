clc;clear;close all;
load("Riskin_Cynopterus_02_th4en5_Es5e4_Ec5e5_Bs38en4_Bc114en4_SFs6en1_SFc8en1_Edge2xpt5xstiff_rho1225en3_SpherePurts_angStdv1Deg_simResults.mat")

meshPoseListCell=permute(meshPoseListCell,[2 3 1]);
meshPoseList=cell2mat(meshPoseListCell);

meshVelListCell=permute(meshVelListCell,[2 3 1]);
meshVelList=cell2mat(meshVelListCell);

numFrames=size(meshPoseList,3);


%% get forces and spring strains
forceTotalMatrix=zeros(numFrames,3);
numSprings=sum(~isnan(L0mat),"all")/2;
inputMatrix=zeros(numFrames,numSprings);
for k=1:numFrames
    meshPoses=meshPoseList(:,:,k); meshVels=meshVelList(:,:,k);

    airForces = getAirForces(meshPoses,meshVels,elemMembership,bladeElemInfo,...
    A,B,C,rho,distribFunc);
    totalAirForce=sum(airForces);
    forceTotalMatrix(k,:)=totalAirForce;

    strains = getSpringStrains(meshPoses,L0mat);
    inputMatrix(k,:)=strains;
end
%% get sensor poses
sensePose=getSpringLocs(refMeshPts,L0mat);
sensePose=sensePose';

%%
x=inputMatrix'; y=forceTotalMatrix';
totalNumSensors=size(x,1);
x=(x-mean(x,2)); y=(y-mean(y,2));
[xDim, numSteps]=size(x); yDim=size(y,1);
SigX=(x*x')/numSteps;
% regularize to avoid singular covariance in x
[Ux,Sx,~]=svd(SigX);
varFloor=totalNumSensors*eps(max(Sx,[],"all"));
Sx=diag(Sx); Sx(Sx<varFloor)=varFloor; Sx=diag(Sx);
SigY=(y*y')/numSteps;
SigXY=(x*y')/numSteps;

MSE=@(B) MSEfun(B,Ux,Sx,SigY,SigXY);

numPairsToCheck=1e3;
numCounterExamps=0;
for k=1:numPairsToCheck
BList=rand(totalNumSensors,2).*rand(1,2)*1e6;
MSEofMeanB=MSE(mean(BList,2));
MeanOfMSEB=mean([MSE(BList(:,1)) MSE(BList(:,2))]);
if MSEofMeanB>MeanOfMSEB
    numCounterExamps=numCounterExamps+1
end

end

if numCounterExamps==0
    disp(['out of ' num2str(numPairsToCheck) ' argument pairs checked, no convexity counter examples found'])
end


%% functions

function [MSE]=MSEfun(B,Ux,Sx,SigY,SigXY)
xDim=length(B);
Ix=eye(xDim);
Bmat=diag(B);
SxInv=diag(1./diag(Sx));

woodburyInv=(SxInv+Ux'*Bmat*Ux)\Ix;
WBinv=SigXY'*(Ix-Bmat*Ux*woodburyInv*Ux');
W=WBinv*Bmat;

MSE=trace(W*Ux*Sx*Ux'*W')-2*trace(W*SigXY)+trace(SigY)+trace(W*WBinv');
end