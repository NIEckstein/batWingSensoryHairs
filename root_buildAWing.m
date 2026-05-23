clc;clear;close all;

% set up a list of the skeletal marker names
skelMrks{1}='shdR';
skelMrks{2}='elbR';
skelMrks{3}='wstR';
skelMrks{4}='t3R';
skelMrks{5}='t4R';
skelMrks{6}='t5R';
skelMrks{7}='mcp3R';
skelMrks{8}='mcp4R';
skelMrks{9}='mcp5R';
skelMrks{10}='ip3R';
skelMrks{11}='ip4R';
skelMrks{12}='ip5R';
skelMrks{13}='hipR';
skelMrks{14}='kneeR';
skelMrks{15}='ankR';
% skelMrks{16}='str'; % ignoring these because theyre inside the body
% skelMrks{17}='lmb'; % will just draw a line from the hip to the shoulder
% to define the medial patagium edge
numSkelMrks=15;
% set up base skeletal connectivity matrix 
skelCmat=false(numSkelMrks);
skelCmat(1,2)=true;% humerus
skelCmat(2,3)=true;% forearm
skelCmat(3,9)=true;% 5th mc
skelCmat(9,12)=true;% 5th prox phalange
skelCmat(12,6)=true;% 5th dist phalange
skelCmat(3,8)=true;% 4th mc
skelCmat(8,11)=true;% 4th prox phalange
skelCmat(11,5)=true;% 4th dist phalange
skelCmat(3,7)=true;% 3rd mc
skelCmat(7,10)=true;% 3rd prox phalange
skelCmat(10,4)=true;% 3rd dist phalange
skelCmat(13,14)=true;% femur
skelCmat(14,15)=true;% foreleg
skelCmat(1,13)=true;% fake bone connecting hip to shoulder
skelCmat=skelCmat|skelCmat';

% set up patagia info: struct with lists of joint marker indices representing vertices of
% the boundary of each patagium, with the first and last indices corresponding
% end points of line segments that form free edges (edges of patagia not fixed to bone)
patInfo.proPat.borderInds=[1;2;3];
patInfo.distDactPat.borderInds=[5;11;8;3;7;10;4];
patInfo.proxDactPat.borderInds=[6;12;9;3;8;11;5];
patInfo.plagiPat.borderInds=[15;14;13;1;2;3;9;12;6];


%% load in data
datFileName='BatData_Riskin_Cynopterus_02';
ogData=load(['batData\' datFileName '.mat']);
numFrames=size(ogData.data,1);
mrkHist=reshape(ogData.data,numFrames,3,[]);
% downSample and trim the data 
downSampNum=3;
startFrame=1;
mrkHist=downsample(mrkHist(startFrame:end,:,:),downSampNum);
numFrames=size(mrkHist,1);
tList=downsample(ogData.t(startFrame:end),downSampNum);
tList=tList-tList(1);

% select out the skeletal markers
skelInds=[]; skelMrksPresent=[];
for skelMrk=1:numSkelMrks
    ind=ogData.mkr.(skelMrks{skelMrk});
    if ~isnan(ind)
        skelInds=[skelInds ind];
        skelMrksPresent=[skelMrksPresent skelMrks(skelMrk)];
    end
end
numPoints=length(skelInds);
mrkHist=mrkHist(:,:,skelInds);

% flip to get it in the stated orientation
mrkHist(:,2,:)=-mrkHist(:,2,:);

% add in translation of the free stream
freeStreamVel=ogData.U;
freeStreamFrameTrans=-tList'*freeStreamVel;
mrkHist(:,1,:)=mrkHist(:,1,:)-freeStreamFrameTrans;

%% get pairwise distance histories for each pair of skeletal markers
% will use these to find the bone lengths and the upper bounds on distances
% between joint pairs that don't share a bone for the reference
% configuration

distHist=zeros(numPoints*(numPoints-1)/2,numFrames);
for frame=1:numFrames
    distHist(:,frame)=pdist(permute(mrkHist(frame,:,:),[3 2 1]));
end
boneInds=squareform(skelCmat);
nonBoneInds=~boneInds;
maxNonBoneDists=max(distHist(nonBoneInds,:),[],2);
minNonBoneDists=min(distHist(nonBoneInds,:),[],2);
boneLengths=mean(distHist(boneInds,:),2);
boneLengthStdvs=std(distHist(boneInds,:),[],2);
fracStdv=boneLengthStdvs./boneLengths;% just to check that the bones are roughly rigid

%% get initial guess for the reference config
% find the frame where the 1st 2 PCs capture the most
% variance, then project down onto those 2 PCs for the initial guess
max2Dvar=0;
for frame=1:numFrames
    [Comps,Scores,Vars]=pca(permute(mrkHist(frame,:,:),[3 2 1]));
    var2D=sum(Vars(1:2));
    if var2D>max2Dvar
        refPts0=Scores(:,1:2);
        max2Dvar=var2D;
    end
end


%% find the reference config
% maximize the total variance, while keeping bone lengths correct and not
% making any pairwise distance between joints larger than some factor times
% its max in-flight value
refNonBoneDistMult=1.1;% controls how splayed the reference wing is;
options = optimoptions('fmincon','MaxFunctionEvaluations',1e8,...
    'MaxIterations',1e5,'EnableFeasibilityMode',true);
refMrkPts=fmincon(@negVarFun,refPts0,[],[],[],[],[],[],@(pts)conFun(pts,boneLengths,refNonBoneDistMult*maxNonBoneDists,minNonBoneDists,boneInds),options);

% visualize the reference points and flip things if necessary
plotARefConfig(refMrkPts,skelCmat);
flipSpan=1;
flipCord=0;
% flipSpan=input('wing tip pointing left? enter 1 for yes, 0 for no: ');
% flipCord=input('leading edge pointing down? enter 1 for yes, 0 for no: ');
if flipSpan||flipCord
    if flipSpan
        refMrkPts(:,1)=-refMrkPts(:,1);
    end
    if flipCord
        refMrkPts(:,2)=-refMrkPts(:,2);
    end
    % plot the corrected frame
    close(gcf)
    plotARefConfig(refMrkPts,skelCmat,'k');
end
hold on;



%% define the base mesh
% set the minimum and maximum allowable spring lengths in the reference
% config
sLengthMax=.008; sLengthMin=.006;
% use matlabs mesh generator to make a triangular mesh adhering to these
% edge length limits
patNames=fieldnames(patInfo);
geomOrder="linear";
refMeshPts=[];elemMembership=[];numMeshPts=0;
for patNum=1:length(patNames)
    borderInds=patInfo.(patNames{patNum}).borderInds;
    borderPts=refMrkPts(borderInds,:);
    pgon=polyshape(borderPts(:,1),borderPts(:,2));
    Tr=triangulation(pgon);
    gm = fegeometry(Tr);
    gm = generateMesh(gm,"Hmax",sLengthMax,"Hmin",sLengthMin,"GeometricOrder",geomOrder);
    meshPtsForThisPat=gm.Mesh.Nodes;
    refMeshPts=[refMeshPts meshPtsForThisPat];
    elemMembership=[elemMembership gm.Mesh.Elements+numMeshPts];
    numMeshPts = size(refMeshPts, 2);
end

refMeshPts=refMeshPts';
elemMembership=elemMembership';
% use elem membership to find connected mesh pts
meshCmat=false(numMeshPts);
for elem=1:size(elemMembership,1)
    conPairs=nchoosek(elemMembership(elem,:),2);
    conPairs=[conPairs;fliplr(conPairs)];% Ensure symmetry in connectivity
    linInds=sub2ind([numMeshPts numMeshPts],conPairs(:,1), conPairs(:,2));
    meshCmat(linInds) = true;
    
end
meshCmatWreps=meshCmat;
% delete duplicate mesh points along the skeleton
tol=.5*sLengthMin/max(abs(refMeshPts(:)));
[refMeshPts,keptInds,IC]= uniquetol(refMeshPts,tol,'ByRows',true);
numMeshPts=size(refMeshPts,1);
% give duplicate points eachother's connections before pruning the cmat
for k=1:length(IC)
    insts=IC(k)==IC;
    numInsts=sum(insts);
    if numInsts>1
        sharedConVec=any(meshCmat(:,insts),2);
        sharedConBlock=repmat(sharedConVec,[1 numInsts]);
        meshCmat(:,insts)=sharedConBlock;
        meshCmat(insts,:)=sharedConBlock';
    end
end
% now prune the meshCmat
meshCmat=meshCmat(keptInds,keptInds);
scatter(refMeshPts(:,1),refMeshPts(:,2));
axis equal
% update the element membership matrix to reflect the new set of
% non-duplicated masses
elemMembership=IC(elemMembership);

%% find edge springs
% first get edge masses
totalWingPolyVerts=refMrkPts([1 3 7 10 4 5 6 15 14 13],:);
[in,on]=inpolygon(refMeshPts(:,1),refMeshPts(:,2),totalWingPolyVerts(:,1),totalWingPolyVerts(:,2));
edgeMassInds=(~in)|on;
% get a cMat for only edge masses
edgeMeshCmat=edgeMassInds & edgeMassInds';

%% define the bone meshes 
% just indices for the frame points and how to construct them from the joint markers
% we'll write functions for their motion later
% also, i'm calling the line from the hip to the shoulder a bone
convThresh=1e-2;
[r, c]=find(tril(skelCmat));
boneMrkIndPairs=[r c];
numBones=length(boneMrkIndPairs);
boneMeshIndsCell=cell(numBones,1);boneWsCell=cell(numBones,1);
for boneNum=1:numBones
    boneEnds = refMrkPts(boneMrkIndPairs(boneNum,:),:);
    w=refMeshPts/boneEnds;
    wSumErr=abs(sum(w,2)-1);
    boneMeshIndsCell{boneNum}=find(~any(w<-convThresh,2)&~any(w>1+convThresh,2)&wSumErr<convThresh);
    w=w(boneMeshIndsCell{boneNum},:);
    % find nearest true conv combo mat to w
    w=findNearestConvComboMat(w);
    boneWsCell{boneNum}=w;
    % use it to clean up the positions of the current bone mesh points
    refMeshPts(boneMeshIndsCell{boneNum},:)=w*boneEnds;
end

% make sure each bone point is only in one bone
bonePairs=nchoosek((1:numBones),2);
for k=1:size(bonePairs,1)
    bonePair=bonePairs(k,:);
    compMat=boneMeshIndsCell{bonePair(1)}==boneMeshIndsCell{bonePair(2)}';
    if any(compMat(:))
        [indsToCutFromBone1,~]=find(compMat);
        boneMeshIndsCell{bonePair(1)}(indsToCutFromBone1)=[];
        boneWsCell{bonePair(1)}(indsToCutFromBone1,:)=[];
        disp('fixed a shared pt across bones')
    end
end
% cut connections between bone mesh points
boneBinInds=any((1:numMeshPts)==cell2mat(boneMeshIndsCell));
boneConnections=boneBinInds&boneBinInds';
meshCmat(boneConnections)=false;
% plot the connections to make sure everything worked
plotARefConfig(refMeshPts,meshCmat,'c');
% plot the bones
hold on;
for boneNum=1:numBones
    bonePts=refMeshPts(boneMeshIndsCell{boneNum},:);
    scatter(bonePts(:,1),bonePts(:,2),'filled');
    hold on;
end
axis equal
% fix the edge spring cmat so it doesn't include connections between
% non-adjacent edge masses
edgeMeshCmat=edgeMeshCmat & meshCmat;

%% set fractional masses based on summed area of the triangles they are involved with 
% (only for masses not in bones, bone masses are nan)
% we'll multiply these by the total patagium mass to get actual masses in
% the sim code
meshFracMasses=nan(numMeshPts,1);
allBoneMeshInds=cell2mat(boneMeshIndsCell);
allPatMeshInds=find(~boneBinInds);
xMeshCords=refMeshPts(:,1); yMeshCords=refMeshPts(:,2);
totalPatArea=0;
for ptNum=1:numMeshPts
    elemsWithThisPt=elemMembership(any(elemMembership==ptNum,2),:);
    x=xMeshCords(elemsWithThisPt);
    y=yMeshCords(elemsWithThisPt);
    if size(elemsWithThisPt,1)==1
        x=x'; y=y';
    end
    AList=abs(sum(x.*[y(:,2)-y(:,3) y(:,3)-y(:,1) y(:,1)-y(:,2)],2));
    patchArea=sum(AList);
    totalPatArea=totalPatArea+patchArea/3;
    if ~any(allBoneMeshInds==ptNum)
       meshFracMasses(ptNum)=patchArea;
    end
end
meshFracMasses=meshFracMasses/sum(meshFracMasses,'omitmissing');
%% define blades 
% set spanwise direction to go from hip to tip in the ref config
refSpanDir=refMrkPts(4,:)-refMrkPts(13,:);
refSpanDir=refSpanDir/norm(refSpanDir);
% rotate and shift reference stuff so hip is at o and axes are span/cord
refCordDir=[-refSpanDir(2) refSpanDir(1)];
rot=[refSpanDir;refCordDir];
refMrkPts=refMrkPts*rot'; 
refMeshPts=refMeshPts*rot';
refMeshPts=refMeshPts-refMrkPts(13,:);
refMrkPts=refMrkPts-refMrkPts(13,:);
% get element centroids  
iList=elemMembership(:,1); 
jList=elemMembership(:,2); 
kList=elemMembership(:,3);
iPts=refMeshPts(iList,:); 
jPts=refMeshPts(jList,:);
kPts=refMeshPts(kList,:);
centers=(iPts+jPts+kPts)/3;
cx=centers(:,1);
% set boundaries between blades
numBlades=10; % from fan & breuer 2021
tipX=refMrkPts(4,1);
xCutCoords=linspace(0,tipX,numBlades+1);
x=refMeshPts(:,1); y=refMeshPts(:,2);
yGE0=y>=-1e-5; yLE0=y<=1e-5;
%allocate elements to blades
for bladeNum=1:numBlades

    xMed=xCutCoords(bladeNum);
    xLat=xCutCoords(bladeNum+1);
    elemInds=find((cx>xMed)&(cx<xLat));
    bladeElemInfo(bladeNum).elemInds=elemInds;

    elems=elemMembership(elemInds,:);
    ptInds=unique(elems(:));
    ptBinInds=false(numMeshPts,1); ptBinInds(ptInds)=true;

    LedgeInds=find(ptBinInds&yGE0&edgeMassInds);
    [~,idxM]=min(x(LedgeInds)); [~,idxL]=max(x(LedgeInds));
    MLind=LedgeInds(idxM); LLind=LedgeInds(idxL);
    

    TedgeInds=find(ptBinInds&yLE0&edgeMassInds);
    [~,idxM]=min(x(TedgeInds)); [~,idxL]=max(x(TedgeInds));
    MTind=TedgeInds(idxM); LTind=TedgeInds(idxL);

    if MTind==MLind
        [~,idx]=min(y(TedgeInds));
        MTind=TedgeInds(idx);
    end
    medInds=[MLind MTind];

    % if LTind~=LLind
    %     latInds=[LLind LTind];
    % else
    %     latInds=LLind;
    % end
    latInds=[LLind LTind];
    
    bladeElemInfo(bladeNum).medInds=medInds;
    bladeElemInfo(bladeNum).latInds=latInds;
end



%% define spring geometric properties 
% (we set material properties in the sim code)

% set unpretensed reference spring lengths (we'll pretense in the sim code)
refSpringLengthMat=squareform(pdist(refMeshPts));% this will get written over later, but in this section, we will use it as is
refSpringLengthMat(~meshCmat)=nan;

% % find spring directions relative to their stiffness principle axes
springDirX=squareform(pdist(refMeshPts(:,1),@minus))./refSpringLengthMat; 
springDirY=squareform(pdist(refMeshPts(:,2),@minus))./refSpringLengthMat;

springDirSpanMat=nan(numMeshPts); springDirCordMat=nan(numMeshPts);
springPatMembershipMat=nan(numMeshPts);
colorsForCheck=['r' 'b' 'g' 'c'];
figure()
for patNum=1:length(patNames)
    % get border points for this pat
    borderInds=patInfo.(patNames{patNum}).borderInds;
    borderPts=refMrkPts(borderInds,:);
    % get span and cord directions for this pat 
    spanDir=borderPts(1,:)-borderPts(end,:);% spanwise direction for each potagium based on the free membrain edge
    spanDir=spanDir/norm(spanDir);
    cordDir=[-spanDir(2) spanDir(1)];% cordwise is just perpendicular to spanwise
    % find springs in this pat
    in=inpolygon(refMeshPts(:,1),refMeshPts(:,2),borderPts(:,1),borderPts(:,2));
    patCmat=(in&in')&meshCmat;
    springPatMembershipMat(patCmat)=patNum;
    % get span-wise spring direction components
    springDirSpanMat(patCmat)=springDirX(patCmat)*spanDir(1)+springDirY(patCmat)*spanDir(2);
    % get cord-wise spring direction components
    springDirCordMat(patCmat)=springDirX(patCmat)*cordDir(1)+springDirY(patCmat)*cordDir(2);
    plotARefConfig_noNewFig(refMeshPts,patCmat,colorsForCheck(patNum))
    hold on;
end
hold off;
%% get purturbed marker trajectories
% % build joint structure to guide pertubation process
joints(1).root=1;% shoulder to arm
joints(1).children=[2 3 9 12 6 8 11 5 7 10 4];
joints(2).root=1;% shoulder to body
joints(2).children=[13 14 15];
joints(3).root=2;% elbow
joints(3).children=[3 9 12 6 8 11 5 7 10 4];
joints(4).root=13;% hip
joints(4).children=[14 15];
joints(5).root=3;% wrist to 5th dig
joints(5).children=[9 12 6];
joints(6).root=3;% wrist to 4th dig
joints(6).children=[8 11 5];
joints(7).root=3;% wrist to 3rd dig
joints(7).children=[7 10 4];
joints(8).root=14;% knee
joints(8).children=15;
joints(9).root=9;% mcp5
joints(9).children=[12 6];
joints(10).root=8;% mcp4
joints(10).children=[11 5];
joints(11).root=7;% mcp3
joints(11).children=[10 4];
joints(12).root=12;% icp5
joints(12).children=6;
joints(13).root=11;% icp4
joints(13).children=5;
joints(14).root=10;% icp3
joints(14).children=4;

numJoints=14;
% % get pruturbed marker trajectories
numPurts=7;% how many purturbed trajectories
numKnots=6;% how many Knots in the angular purturbation spline, including end points
tDistLim=tList(end)/(numKnots-1)/3;% two knots can't be less than this close in time, prevents extreme accs
angStdv=1*pi/180;% roughly the standard deviation of the angle we rotate at each joint to cause the purturbation
mrkHistPurts=zeros([size(mrkHist) numPurts]);
for purtNum=1:numPurts

    % % get smooth sequences of  unit quaternions for each joint
    minTdist=0;
    while minTdist<tDistLim
        knots=[0 tList(end)];
        knots=sort([knots rand(1, numKnots-2)*tList(end)]);
        minTdist=min(diff(knots));
    end

    knotDirs=randn(numJoints,3,numKnots)*(angStdv/2/sqrt(3)).*sqrt(3*[4/5 1/10 1/10]);% comment/uncomment ".*sqrt(3*[4/5 1/10 1/10])" at the end of this line to toggle between spherically symmetric joint angle purturbations and purturbations that are biased to rotate around the free stream velocity axis
    dirs=spline(knots,knotDirs,tList);
    quats=[ones(numJoints,1,length(tList)) dirs];
    quats=quats./vecnorm(quats,2,2);

    % % do the purturbation 
    for tStep=1:length(tList)
        mrksNow=mrkHist(tStep,:,:);
        quatsNow=quats(:,:,tStep);
        for jNum=1:numJoints
            rot=quat2rotm(quatsNow(jNum,:));
            rootInd=joints(jNum).root;
            childInds=joints(jNum).children;

            root=mrksNow(:,:,rootInd);
            children=mrksNow(:,:,childInds);

            children=children-root;
            children=pagemtimes(children,rot);
            children=children+root;
            mrksNow(:,:,childInds)=children;
        end
        mrkHistPurts(tStep,:,:,purtNum)=mrksNow;
    end

end
% stick the original on the front of the purturbations
mrkHistPurts=cat(4,mrkHist,mrkHistPurts);
% % animate to check that nothing's crazy
% reorient for plotting 
mrkHistPurts=permute(mrkHistPurts,[3 2 1 4]);
colors=hsv(numPurts+1);
[indsFrom,indsTo]=find(tril(skelCmat));
buffer=.05;
xLims=[min(mrkHistPurts(:,1,:,:),[],"all")-buffer max(mrkHistPurts(:,1,:,:),[],"all")+buffer];
yLims=[min(mrkHistPurts(:,2,:,:),[],"all")-buffer max(mrkHistPurts(:,2,:,:),[],"all")+buffer];
zLims=[min(mrkHistPurts(:,3,:,:),[],"all")-buffer max(mrkHistPurts(:,3,:,:),[],"all")+buffer];
figure()
for purtNum=1:numPurts+1
    mrksNow=mrkHistPurts(:,:,1,purtNum);
    x=mrksNow(:,1)'; y=mrksNow(:,2)'; z=mrksNow(:,3)';
    xPlot=[x(indsFrom);x(indsTo)];
    yPlot=[y(indsFrom);y(indsTo)];
    zPlot=[z(indsFrom);z(indsTo)];
    H(:,purtNum)=plot3(xPlot,yPlot,zPlot,'Color',colors(purtNum,:));
    hold on;
end
grid on;
axis equal;
xlim(xLims); ylim(yLims); zlim(zLims);
for tStep=1:length(tList)
    for purtNum=1:numPurts+1
        mrksNow=mrkHistPurts(:,:,tStep,purtNum);
        x=mrksNow(:,1)'; y=mrksNow(:,2)'; z=mrksNow(:,3)';
        xPlot=[x(indsFrom);x(indsTo)];
        yPlot=[y(indsFrom);y(indsTo)];
        zPlot=[z(indsFrom);z(indsTo)];
        for seg=1:size(H,1)
            H(seg,purtNum).XData=xPlot(:,seg);
            H(seg,purtNum).YData=yPlot(:,seg);
            H(seg,purtNum).ZData=zPlot(:,seg);
        end
    end
    drawnow;
    pause(.05)
end
% put marker history array back in the starting orientation
mrkHistPurts=permute(mrkHistPurts,[3 2 1 4]);

%% make function handels that take in time and output the joint marker positions and velocities
mrkPoseFunList=cell(numPurts+1,1);
mrkVelFunList=mrkPoseFunList;

for purtNum=1:numPurts+1
    posePPcell=cell(numSkelMrks,3);
    velPPcell=posePPcell;
    mrkHistNow=mrkHistPurts(:,:,:,purtNum);
    for dim=1:3
        for mrk=1:numSkelMrks
            coordList=mrkHistNow(:,dim,mrk);
            posePP=spline(tList,coordList);
            posePPcell{mrk,dim}=posePP;
            velPP=posePP;
            velPP.order=velPP.order-1;% derivative drops order
            velPP.coefs=velPP.coefs(:,1:3).*[3 2 1];
            velPPcell{mrk,dim}=velPP;
        end
    end

    mrkPoseFun=@(t) cellfun(@(pp)ppval(pp,t),posePPcell);
    mrkPoseFunList{purtNum}=mrkPoseFun;

    mrkVelFun=@(t) cellfun(@(pp)ppval(pp,t),velPPcell);
    mrkVelFunList{purtNum}=mrkVelFun;
end
%% set up a rough initial batwing state for each purturbation
triangles=delaunay(refMrkPts);
numTriangles=size(triangles,1);
for purtNum=1:numPurts+1
    mrkPoses0=mrkPoseFunList{purtNum}(0);
    mrkVels0=mrkVelFunList{purtNum}(0);
    meshPoses0=zeros(numMeshPts,3);
    meshVels0=meshPoses0;
    for tri=1:numTriangles
        vertInds=triangles(tri,:);
        verts=refMrkPts(vertInds,:);
        indsInTri = inpolygon(refMeshPts(:,1),refMeshPts(:,2),verts(:,1),verts(:,2));
        refMeshPtsInTri=refMeshPts(indsInTri,:);
        numIndsInTri=sum(indsInTri);
        % get conv combo weights
        w=[refMeshPtsInTri ones(numIndsInTri,1)]/[verts ones(3,1)];% ones in the last columns constrain the combo weights to sum to 1
        % get initial poses
        meshPoses0(indsInTri,:)=w*mrkPoses0(vertInds,:);
        % get velocities
        meshVels0(indsInTri,:)=w*mrkVels0(vertInds,:);
    end
    meshPoses0List{purtNum}=meshPoses0;
    meshVels0List{purtNum}=meshVels0;
end
% plot the last one to check
figure()
colors=hsv(numBlades);
for bladeNum=1:numBlades
    elems=elemMembership(bladeElemInfo(bladeNum).elemInds,:);
    patch('Faces',elems,'Vertices',meshPoses0,'FaceColor',colors(bladeNum,:));
    hold on;

    medInds=bladeElemInfo(bladeNum).medInds;
    latInds=bladeElemInfo(bladeNum).latInds;
    medPts=meshPoses0(medInds,:); latPts=meshPoses0(latInds,:);
    scatter3(medPts(:,1), medPts(:,2),medPts(:,3),100,'ro')
    scatter3(latPts(:,1), latPts(:,2),latPts(:,3),200,'go')

end
axis equal;
grid on;
hold off;

%% use initial mesh pose of the unpurturbed version to get our ref spring lengths

refSpringLengthMat=squareform(pdist(meshPoses0List{1}));
refSpringLengthMat(~meshCmat)=nan;


% save out the necessary things to feed into the simulation code/ to remember the key parameters used by this code
outputFileName=[datFileName '_ds3_sRng6en3_8en3_flapBiasedPurts_angStdv1Deg_setUpForSim.mat'];
save(outputFileName,'mrkVelFunList','mrkPoseFunList','boneMeshIndsCell','boneMrkIndPairs',...
    'refSpringLengthMat','refMeshPts','meshFracMasses','boneWsCell','totalPatArea',...
    'elemMembership','allBoneMeshInds','refNonBoneDistMult','sLengthMin',...
    'sLengthMax','geomOrder','convThresh','refMrkPts','tList','meshVels0List',...
    'meshPoses0List','allPatMeshInds','edgeMeshCmat','springDirCordMat',...
    'springDirSpanMat','springPatMembershipMat','bladeElemInfo','freeStreamVel');

%% save out cmats needed for making nice plots
% save('cmats.mat',"meshCmat","skelCmat");

%% functions
function [cineq,ceq]=conFun(pts,refBoneLengths,maxNonBoneDists,minNonBoneDists,boneInds)

dists=pdist(pts);
boneLengths=dists(boneInds);
nonBoneDists=dists(~boneInds);

cineq=[nonBoneDists-maxNonBoneDists';minNonBoneDists'-nonBoneDists];
ceq=boneLengths-refBoneLengths';


end

function neg2Dvar=negVarFun(pts)

neg2Dvar=-trace(cov(pts));

end