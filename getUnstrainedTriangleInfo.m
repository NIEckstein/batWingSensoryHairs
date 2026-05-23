function [jkPoseList] = getUnstrainedTriangleInfo(refMeshPoses,elemMembership,meanSF)
% finds coordinates of the j and k vertices of the reference traingle when 
% the i point is at the origin

iList=elemMembership(:,1);
jList=elemMembership(:,2);
kList=elemMembership(:,3);

refMeshPoses=meanSF*refMeshPoses;

iPoses=refMeshPoses(iList,:); 
jPoses=refMeshPoses(jList,:);jPoses=jPoses-iPoses; 
kPoses=refMeshPoses(kList,:);kPoses=kPoses-iPoses; 

jkPoseList=cat(3,jPoses,kPoses);
jkPoseList=permute(jkPoseList,[3 2 1]);




end