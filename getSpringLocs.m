function Locs = getSpringLocs(refMeshPoses,L0mat)

% setup some things
numMasses=size(refMeshPoses,1);
[iList,jList]=find(tril(~isnan(L0mat)));


% get lists of strains
iPoses=refMeshPoses(iList,:); jPoses=refMeshPoses(jList,:);
Locs=(iPoses+jPoses)/2;


end