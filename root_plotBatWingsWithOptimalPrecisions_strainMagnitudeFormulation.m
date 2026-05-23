clc; clear; close all;

allFileNamesToPlot={'Riskin_Cynopterus_02_th4en5_Es5e4_Ec5e5_Bs38en4_Bc114en4_SFs6en1_SFc8en1_Edge2xpt5xstiff_rho1225en3_SpherePurts_angStdv1Deg_optResults_3DrsqrFrac99en2.mat';...
    'Riskin_Cynopterus_02_th4en5_Es5e4_Ec5e5_Bs38en4_Bc114en4_SFs6en1_SFc8en1_Edge2xpt5xstiff_rho1225en3_SpherePurts_angStdv1Deg_optResults_3DrsqrFrac995en3_2.mat';...
    'Riskin_Cynopterus_02_th4en5_Es5e4_Ec5e5_Bs38en4_Bc114en4_SFs6en1_SFc8en1_Edge2xpt5xstiff_rho1225en3_SpherePurts_angStdv1Deg_optResults_3DrsqrFrac995en3.mat';...
    'Riskin_Cynopterus_02_th4en5_Es5e4_Ec5e5_Bs38en4_Bc114en4_SFs6en1_SFc8en1_Edge2xpt5xstiff_rho1225en3_flapBiasedPurts_angStdv1Deg_optResults_3DrsqrFrac995en3.mat';...
    'Riskin_Cynopterus_02_th4en5_Es25e3_Ec25e4_Bs38en4_Bc114en4_SFs6en1_SFc8en1_Edge2xpt5xstiff_rho1225en3_SpherePurts_angStdv1Deg_optResults_3DrsqrFrac995en3.mat'};
numFiles=length(allFileNamesToPlot);
BoptList=[]; RsqrBoundList=[];
for k=1:numFiles
    load(allFileNamesToPlot{k});
    BoptList=[BoptList Bopt];
    RsqrBoundList=[RsqrBoundList RsqrBound];
end

% [RsqrBoundList,idx]=sort(RsqrBoundList);
% BoptList=BoptList(:,idx);
maxes=max(BoptList);
maxAll=max(maxes);
BoptList=[BoptList;zeros(1,numFiles);maxAll*ones(1,numFiles)];
mrkX=refMrkPts(:,1); mrkY=refMrkPts(:,2);
boneLinesX=mrkX(boneMrkIndPairs'); boneLinesY=mrkY(boneMrkIndPairs');
refMeshPts=[refMeshPts;nan(3,2)];
elemMembership=[elemMembership;[1:3;1:3]+numMeshPts];
cMap=interp1([1;1000],[1 1 1;1 0 0],(1:1000)');
for k=1:numFiles
    %plot the wing with frame filled
    figure();
    colormap(cMap);
    %plot the precisions between the gridPts
    P=patch('Faces',elemMembership,'Vertices',refMeshPts);
    P.FaceColor='flat';
    P.CData=BoptList(:,k);
    P.EdgeColor='none';
    %P.FaceAlpha='flat';
    %P.AlphaDataMapping="none";
    %P.FaceVertexAlphaData=ones(size(BoptList(:,k)));
    colorbar('Limits',[0 maxes(k)]);
    hold on;
    scatter(refMeshPts(allBoneMeshInds,1),refMeshPts(allBoneMeshInds,2),'b','filled');
    plot(boneLinesX,boneLinesY,'b')
    axis equal
    title(['R^2 constraint=' num2str(RsqrBoundList(k)) ', max precis=' num2str(maxes(k))]);

    
    hold off;
end
 