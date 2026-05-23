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


% check MSE with almost no noise to get a baseline and select the MSE bound
RsqrFracBound=.99;
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

% %% check that the best mse achievable when noise corresponding to Bopt is
% % added matches the MSE bound
% numReps=100;
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
precisPlot=sqrt(Bopt/max(Bopt));
figure()
P=scatter(sensePose(1,:),sensePose(2,:),'filled');
P.MarkerFaceAlpha='flat';
P.AlphaData=precisPlot;
%P.AlphaDataMapping="direct";
hold on;
scatter(sensePose(1,:),sensePose(2,:))
axis equal

%%

save('Riskin_Cynopterus_02_th4en5_Es5e4_Ec5e5_Bs38en4_Bc114en4_SFs6en1_SFc8en1_Edge2xpt5xstiff_rho1225en3_SpherePurts_angStdv1Deg_optResults_3DrsqrFrac99en2_uniAx.mat')






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

