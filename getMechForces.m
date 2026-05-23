function Forces = getMechForces(meshPoses,meshVels,Kmat,Bmat,L0mat)

% setup some things needed for both elastic and damping forces
numMasses=size(meshPoses,1);
cMat=tril(~isnan(Kmat));
springlessMasses=find(~any(cMat,2));
[iList,jList]=find(cMat);
linInds=sub2ind([numMasses numMasses],iList,jList);

% get lists of all elastic forces with a sepereate entry per spring/mass
% pair
iPoses=meshPoses(iList,:); jPoses=meshPoses(jList,:);
L0s=L0mat(linInds);
Ks=Kmat(linInds);
diffs=jPoses-iPoses;
lengths=vecnorm(diffs,2,2);
deltas=lengths-L0s;
dirs=diffs./lengths;
iForces=dirs.*Ks.*deltas;
jForces=-iForces;
unCombinedElasticForces=[iForces;jForces];

% get lists of all damping forces with a sepereate entry per spring/mass
% pair
iVels=meshVels(iList,:); jVels=meshVels(jList,:);
Bs=Bmat(linInds);
diffs=jVels-iVels;
mags=sum(diffs.*dirs,2);
iForces=dirs.*Bs.*mags;
jForces=-iForces;
unCombinedDampingForces=[iForces;jForces];

unCombinedForces=unCombinedElasticForces+unCombinedDampingForces;

% combine the forces so we get one summed force per mass
indList=[iList;jList;springlessMasses];
unCombinedForces=[unCombinedForces;zeros(length(springlessMasses),3)];
Forces=splitapply(@(x)sum(x,1),unCombinedForces,indList);


end