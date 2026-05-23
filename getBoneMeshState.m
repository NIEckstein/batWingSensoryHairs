function [boneMeshPoses,boneMeshVels]=getBoneMeshState(t,mrkPoseFun,mrkVelFun,boneMrkIndPairs,boneWsCell)

mrkPoses=mrkPoseFun(t); mrkVels=mrkVelFun(t);
numBones=size(boneMrkIndPairs,1);
numBoneMeshPtsList=cellfun(@(x)size(x,1),boneWsCell);
numBoneMeshPtsTotal=sum(numBoneMeshPtsList);
boneMeshPoses=zeros(numBoneMeshPtsTotal,3);
boneMeshVels=boneMeshPoses;
lastEntry=0;
for boneNum=1:numBones
    mrkInds=boneMrkIndPairs(boneNum,:);

    posePair=mrkPoses(mrkInds,:);
    velPair=mrkVels(mrkInds,:);

    meshPoses=boneWsCell{boneNum}*posePair;
    meshVels=boneWsCell{boneNum}*velPair;

    numMeshNow=numBoneMeshPtsList(boneNum);
    boneMeshPoses(lastEntry+1:lastEntry+numMeshNow,:)=meshPoses;
    boneMeshVels(lastEntry+1:lastEntry+numMeshNow,:)=meshVels;

    lastEntry=lastEntry+numMeshNow;
end

end