function Locs = getSenseLocs(refMeshPoses,elemMembership)

% setup some things

iList=elemMembership(:,1);
jList=elemMembership(:,2);
kList=elemMembership(:,3);


% get lists of strains
iPoses=refMeshPoses(iList,:); 
jPoses=refMeshPoses(jList,:); 
kPoses=refMeshPoses(kList,:); 
Locs=(iPoses+jPoses+kPoses)/3;


end