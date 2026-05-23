clc; clear; close all;
% you cannot compare across strain formulations in this code
%% load in the optimal precisions to compare for strain magnitude formulation
allFileNamesToCompare={'Riskin_Cynopterus_02_th4en5_Es5e4_Ec5e5_Bs38en4_Bc114en4_SFs6en1_SFc8en1_Edge2xpt5xstiff_rho1225en3_SpherePurts_angStdv1Deg_optResults_3DrsqrFrac99en2.mat';...
    'Riskin_Cynopterus_02_th4en5_Es5e4_Ec5e5_Bs38en4_Bc114en4_SFs6en1_SFc8en1_Edge2xpt5xstiff_rho1225en3_SpherePurts_angStdv1Deg_optResults_3DrsqrFrac995en3_2.mat';...
    'Riskin_Cynopterus_02_th4en5_Es5e4_Ec5e5_Bs38en4_Bc114en4_SFs6en1_SFc8en1_Edge2xpt5xstiff_rho1225en3_SpherePurts_angStdv1Deg_optResults_3DrsqrFrac995en3.mat';...
    'Riskin_Cynopterus_02_th4en5_Es5e4_Ec5e5_Bs38en4_Bc114en4_SFs6en1_SFc8en1_Edge2xpt5xstiff_rho1225en3_flapBiasedPurts_angStdv1Deg_optResults_3DrsqrFrac995en3.mat';...
    'Riskin_Cynopterus_02_th4en5_Es25e3_Ec25e4_Bs38en4_Bc114en4_SFs6en1_SFc8en1_Edge2xpt5xstiff_rho1225en3_SpherePurts_angStdv1Deg_optResults_3DrsqrFrac995en3.mat'};
names=["base, low R^2","base, repeat","base","flap-biased perturbations","half stiffness" ];

% %% load in the optimal precisions to compare for uniaxial formulation
% allFileNamesToCompare={'Riskin_Cynopterus_02_th4en5_Es5e4_Ec5e5_Bs38en4_Bc114en4_SFs6en1_SFc8en1_Edge2xpt5xstiff_rho1225en3_SpherePurts_angStdv1Deg_optResults_3DrsqrFrac99en2_uniAx.mat';...
%     'Riskin_Cynopterus_02_th4en5_Es5e4_Ec5e5_Bs38en4_Bc114en4_SFs6en1_SFc8en1_Edge2xpt5xstiff_rho1225en3_SpherePurts_angStdv1Deg_optResults_3DrsqrFrac995en3_uniAx.mat'};
% names=["low R^2","high R^2"];

 
numFiles=length(allFileNamesToCompare);
load(allFileNamesToCompare{1},'Bopt','sensePose','sLengthMin','sLengthMax');
sensePose=sensePose';
BoptList=Bopt;
for k=2:numFiles
    load(allFileNamesToCompare{k},'Bopt');
    BoptList=[BoptList Bopt];
end
clear Bopt;
%% check their pairwise spatial correlation via Lee's L statistic
numPerms=1e5;% number of spatial permutations for getting p-values
neighborRadius=2*(sLengthMax+sLengthMin)/2;% radius that determins neighboring sensors for Lee's L
weightFun=@(V) double(V<neighborRadius);
LMat=zeros(numFiles); pValMat=LMat;
for f1=1:numFiles-1
    for f2=f1+1:numFiles
        BsToCompare=BoptList(:,[f1 f2]);
        [L,pVal]=myLeesL(sensePose,BsToCompare,weightFun,numPerms);
        LMat(f1,f2)=L; LMat(f2,f1)=L;
        pValMat(f1,f2)=pVal; pValMat(f2,f1)=pVal;
    end
end
LMat(logical(eye(numFiles)))=1;


figure()
subplot(1,2,1)
G(1)=heatmap(names,names,LMat);
colormap(G(1),redbluecmap);
G(1).ColorLimits=[-1 1];
title("Lee's L")
subplot(1,2,2)
G(2)=heatmap(names,names,pValMat);
colormap(G(2),flipud(sky));
G(2).ColorLimits=[0 1];
title("one-sided p-Values")