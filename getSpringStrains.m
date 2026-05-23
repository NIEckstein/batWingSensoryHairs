function strains = getSpringStrains(meshPoses,L0mat)

% setup some things
numMasses=size(meshPoses,1);
[iList,jList]=find(tril(~isnan(L0mat)));
linInds=sub2ind([numMasses numMasses],iList,jList);

% get lists of strains
iPoses=meshPoses(iList,:); jPoses=meshPoses(jList,:);
L0s=L0mat(linInds);
diffs=jPoses-iPoses;
lengths=vecnorm(diffs,2,2);
strains=((lengths-L0s)./L0s)';


end