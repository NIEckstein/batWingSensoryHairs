function [pStrains] = getPrincipleStrains(meshPoses,refjkPoseList,elemMembership)

% get jk pose list for all elements
iList=elemMembership(:,1); 
jList=elemMembership(:,2);
kList=elemMembership(:,3);
iPoses=meshPoses(iList,:); 
jPoses=meshPoses(jList,:); jPoses=jPoses-iPoses; 
kPoses=meshPoses(kList,:); kPoses=kPoses-iPoses;
jkPoseList=cat(3,jPoses,kPoses);
jkPoseList=permute(jkPoseList,[3 2 1]);
% find singular values of the transformation that
% goes from the ref pose to the current pose for all elements. 
% These are precisely the principle strains plus 1 
wList=pagelsqminnorm(refjkPoseList,jkPoseList);
svList=pagesvd(wList);
pStrains=svList(1:2,:)'-1;%Selects out first two svs (the third is 0), reshapes to be numElems x 2, and subtracts 1 to get principle strains 



end