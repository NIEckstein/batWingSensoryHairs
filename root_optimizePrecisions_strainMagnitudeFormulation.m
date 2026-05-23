clc;clear;close all;
load("Riskin_Cynopterus_02_th4en5_Es5e4_Ec5e5_Bs38en4_Bc114en4_SFs6en1_SFc8en1_Edge2xpt5xstiff_rho1225en3_SpherePurts_angStdv1Deg_simResults_2.mat")

meshPoseListCell=permute(meshPoseListCell,[2 3 1]);
meshPoseList=cell2mat(meshPoseListCell);

meshVelListCell=permute(meshVelListCell,[2 3 1]);
meshVelList=cell2mat(meshVelListCell);

numFrames=size(meshPoseList,3);

%% get forces and strain measurments
forceTotalMatrix=zeros(numFrames,3);
numElems=size(elemMembership,1);
inputMatrix=zeros(numFrames,numElems);
meanSF=(ssfSpan+ssfCord)/2;
refjkPoseList = getUnstrainedTriangleInfo(meshPoses0List{1},elemMembership,meanSF);


for k=1:numFrames
    meshPoses=meshPoseList(:,:,k); meshVels=meshVelList(:,:,k);

    airForces = getAirForces(meshPoses,meshVels,elemMembership,bladeElemInfo,...
    A,B,C,rho,distribFunc);
    totalAirForce=sum(airForces);
    forceTotalMatrix(k,:)=totalAirForce;
    % get principle strains
    [pStrains]= getPrincipleStrains(meshPoses,refjkPoseList,elemMembership);
    % get total strain measurment (im just doing 2-norm of principle strains, we could handle this more efficiently without ever computing the svd, but better to have the flexibility to do just dilation or distortion measures if we choose)
    meas=vecnorm(pStrains,2,2);
    inputMatrix(k,:)=meas;
end

%% get sensor poses
sensePose=getSenseLocs(refMeshPts,elemMembership);
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


% check MSE with almost no noise to get a baseline and select the MSE bound
RsqrFracBound=.995;
BlowNoise=ones(xDim,1)/varFloor;
MSElowNoise=MSEfun(BlowNoise,Ux,Sx,SigY,SigXY);
RsqrLowNoise=1-MSElowNoise/trace(SigY)
RsqrBound=RsqrLowNoise*RsqrFracBound
MSEbound=(1-RsqrBound)*trace(SigY);
%% get a good initial guess
% find the precision that gets us to the MSE bound if all sensors have it
initCost=@(p)(MSEfun(abs(p)*ones(xDim,1),Ux,Sx,SigY,SigXY)-MSEbound)^2;
allStartPrecis=fminsearch(initCost,1000);

%%
tic;
% run the optimization
% options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',1e6,...
%     'MaxIterations',1e4,'ScaleProblem',true,'UseParallel',true);
options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',1e8,...
    'MaxIterations',1e5,'UseParallel',true,'EnableFeasibilityMode',true);
% options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',1e6,...
%     'MaxIterations',1e4,'UseParallel',true);


Bopt=fmincon(@costFun,allStartPrecis*ones(size(BlowNoise)),[],[],[],[],zeros(xDim,1),[],@(B)MSEcon(B,MSEbound,Ux,Sx,SigY,SigXY),options);
runTime=toc;
% % run it again from a different IC to check convexity
% Bopt2=fmincon(@costFun,1e5*abs(randn(xDim,1)),[],[],[],[],zeros(xDim,1),[],@(B)MSEcon(B,MSEbound,Ux,Sx,SigY,SigXY),options);

%% check that the best mse achievable when noise corresponding to Bopt is
% added matches the MSE bound
% numReps=1;
% noise=randn(xDim,numSteps*numReps).*sqrt(1./Bopt);
% xTest=repmat(x,[1 numReps])+noise;
% yTest=repmat(y,[1 numReps]);
% Wtest=yTest/xTest;
% MSEtest=mean(sum((yTest-Wtest*xTest).^2,1));
% figure()
% plot((Wtest*xTest)');
% hold on
% plot(yTest');
% legend('model','true')

%% plot the precisions
% precisPlot=sqrt(Bopt/max(Bopt));
% figure()
% P=scatter(sensePose(1,:),sensePose(2,:),'filled');
% P.MarkerFaceAlpha='flat';
% P.AlphaData=precisPlot;
% %P.AlphaDataMapping="direct";
% hold on;
% scatter(sensePose(1,:),sensePose(2,:))
% axis equal



save('Riskin_Cynopterus_02_th4en5_Es5e4_Ec5e5_Bs38en4_Bc114en4_SFs6en1_SFc8en1_Edge2xpt5xstiff_rho1225en3_SpherePurts_angStdv1Deg_optResults_3DrsqrFrac995en3_2.mat.mat')





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

function [cineq,ceq]=MSEcon(B,MSEbound,Ux,Sx,SigY,SigXY)

MSE=MSEfun(B,Ux,Sx,SigY,SigXY);

cineq=MSE-MSEbound;
ceq=[];
end

function cost=costFun(B)

cost=sum(abs(B));

end

