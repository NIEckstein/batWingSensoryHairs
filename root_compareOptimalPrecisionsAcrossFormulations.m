clc; clear; close all;

%% load in optimization data (you must compare results from both strain formulations for this code)
% % List the optimal precision files to compare for strain magnitude formulation
allStrainMagFileNamesToCompare={'Riskin_Cynopterus_02_th4en5_Es5e4_Ec5e5_Bs38en4_Bc114en4_SFs6en1_SFc8en1_Edge2xpt5xstiff_rho1225en3_SpherePurts_angStdv1Deg_optResults_3DrsqrFrac995en3.mat';...
    'Riskin_Cynopterus_02_th4en5_Es5e4_Ec5e5_Bs38en4_Bc114en4_SFs6en1_SFc8en1_Edge2xpt5xstiff_rho1225en3_SpherePurts_angStdv1Deg_optResults_3DrsqrFrac99en2.mat';...
    'Riskin_Cynopterus_02_th4en5_Es5e4_Ec5e5_Bs38en4_Bc114en4_SFs6en1_SFc8en1_Edge2xpt5xstiff_rho1225en3_SpherePurts_angStdv1Deg_optResults_3DrsqrFrac995en3_2.mat';...
    'Riskin_Cynopterus_02_th4en5_Es5e4_Ec5e5_Bs38en4_Bc114en4_SFs6en1_SFc8en1_Edge2xpt5xstiff_rho1225en3_flapBiasedPurts_angStdv1Deg_optResults_3DrsqrFrac995en3.mat';...
    'Riskin_Cynopterus_02_th4en5_Es25e3_Ec25e4_Bs38en4_Bc114en4_SFs6en1_SFc8en1_Edge2xpt5xstiff_rho1225en3_SpherePurts_angStdv1Deg_optResults_3DrsqrFrac995en3.mat'};

% % List the optimal precisions files to compare for uniaxial formulation
allUniAxFileNamesToCompare={'Riskin_Cynopterus_02_th4en5_Es5e4_Ec5e5_Bs38en4_Bc114en4_SFs6en1_SFc8en1_Edge2xpt5xstiff_rho1225en3_SpherePurts_angStdv1Deg_optResults_3DrsqrFrac995en3_uniAx.mat';...
    'Riskin_Cynopterus_02_th4en5_Es5e4_Ec5e5_Bs38en4_Bc114en4_SFs6en1_SFc8en1_Edge2xpt5xstiff_rho1225en3_SpherePurts_angStdv1Deg_optResults_3DrsqrFrac99en2_uniAx.mat'};
names=["strain magnitude", "strain magnitude, lower R^2","strain magnitude, new pert realization","strain magnitude, flap-biased perturbations","strain magnitude, half stiffness","uniaxial strain","uniaxial strain, lower R^2"];

% % load in the strain magnidude formulation results
numStrainMagFiles=length(allStrainMagFileNamesToCompare);
BoptList=[];
for k=1:numStrainMagFiles
    load(allStrainMagFileNamesToCompare{k},'Bopt','sensePose','sLengthMin','sLengthMax');
    BoptList=[BoptList Bopt];
end
% % load in the uniaxial formulation results
numUniAxFiles=length(allUniAxFileNamesToCompare);
for k=1:numUniAxFiles
    load(allUniAxFileNamesToCompare{k},'Bopt','L0mat','elemMembership');
    Bopt=uniax2mag(Bopt,~isnan(L0mat),elemMembership);
    BoptList=[BoptList Bopt];
end

numFiles=numStrainMagFiles+numUniAxFiles;
sensePose=sensePose';
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

%% Functions

function [BsMag] = uniax2mag(BuniAx,cmat,elemMembership)
% converts optimal strains from uniaxial formulation to strain mag
% formulation by assigning the mean of the precisions on the edges of each
% triangle element to the hypothetical strain magnitude sensor associated
% with that triangle
BuniAxMat=nan(size(cmat));
BuniAxMat(tril(cmat))=BuniAx;
BsMag=zeros(size(elemMembership,1),1);
for elem=1:size(elemMembership,1)
    vertInds=elemMembership(elem,:);
    pairs=nchoosek(vertInds,2);
    pairs=sort(pairs,2,"descend");
    linInds=sub2ind(size(cmat),pairs(:,1),pairs(:,2));
    BsMag(elem)=mean(BuniAxMat(linInds),"omitmissing");
end
% not all strain mag sensors have uniaxial sensors on their edges, so we
% fill these values with the median
BsMag(isnan(BsMag))=median(BsMag,"omitmissing");

end