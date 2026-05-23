function Forces = getAirForces(meshPoses,meshVels,elemMembership,bladeElemInfo,...
    A,B,C,rho,distribFunc)

% get indices of first second and third triangle verticies seperated
iList=elemMembership(:,1); jList=elemMembership(:,2); kList=elemMembership(:,3);

% get triangle areas and centroids
iPoses=meshPoses(iList,:); jPoses=meshPoses(jList,:); kPoses=meshPoses(kList,:);
v1List=jPoses-iPoses; v2List=kPoses-jPoses;
cVecs=cross(v1List,v2List,2);
Ax2List=vecnorm(cVecs,2,2);
AList=Ax2List/2;% areas
centPoses=(iPoses+jPoses+kPoses)/3;

%% loop through the blades to get forces per triangle element
% we'll distribute to masses after the loop with splitapply
numBlades=length(bladeElemInfo);
numTris=length(elemMembership);
triangleForces=zeros(numTris,3);
options = optimoptions('fmincon','display','none',...
    'SpecifyConstraintGradient',true,'SpecifyObjectiveGradient',true);
for bNum=1:numBlades
    % get things for this blade
    medInds=bladeElemInfo(bNum).medInds; latInds=bladeElemInfo(bNum).latInds;
    elemInds=bladeElemInfo(bNum).elemInds;
    cents=centPoses(elemInds,:);
    As=AList(elemInds);
    %% get the triangles in their local coordinate frame, origin at the LE
    % get local xy frame
    medCorners=meshPoses(medInds,:); latCorners=meshPoses(latInds,:);
    xDir=-(diff(medCorners)+diff(latCorners))/2;
    yDir= mean(medCorners)-mean(latCorners); 
    % Cost=@(Vars) orthoCost(Vars,[xDir yDir]');
    % Vars=fmincon(Cost,[xDir yDir]',[],[],[],[],[],[],@orthoCon,options);
    %xDir=Vars(1:3)'; yDir=Vars(4:6)';
    xDir=xDir/norm(xDir);
    yDir=yDir*(eye(3)-xDir'*xDir);% project yDir into space orthogonal to xDir (we're trusting xDir more than yDir since it's longer)
    yDir=yDir/norm(yDir);
    zDir=cross(xDir,yDir);
    rotMat=[xDir;yDir;zDir]';
    % put triangle centers in local frame
    centsLocal=cents*rotMat;
    [~,oInd]=max(centsLocal(:,1));
    centsLocal=centsLocal-centsLocal(oInd,:);
    %% get cord length and spanwise thickness of this blade
    spanThick=range(centsLocal(:,2));
    cordL=range(centsLocal(:,1));

    %% get zero lift AoA
    % put x and z in units of cordL, and have x increasing from leading to
    % trailing
    xRel=-centsLocal(:,1)/cordL;
    zRel=centsLocal(:,3)/cordL;
    % fit a quadratic bezier with 1st/3rd control points clamped to
    % leading/trailing edges
    [maxCam,maxCamLoc] = fitBezierByLM(xRel,zRel,10);% third argument is how many Levenberg-Marquardt steps to take
    alph0=-atan2(maxCam,1-maxCamLoc);% zero-lift AoA, Supercool/Jones method https://www.supercoolprops.com/home/articles/camber_and_zero_lift_angle.html

    %% get effective AoA
    % get quarter cord velocity in local blade coords
    medCornerVels=meshVels(medInds,:); latCornerVels=meshVels(latInds,:);
    leadVel=(medCornerVels(1,:)+latCornerVels(1,:))/2;
    trailVel=(medCornerVels(2,:)+latCornerVels(2,:))/2;
    qCordVel= .75*leadVel+.25*trailVel;
    % put it in ,ocal coords to get the AoA
    qCordVelLocal=qCordVel*rotMat;
    alph=-atan2(qCordVelLocal(3),qCordVelLocal(1));% effective AoA

    %%  get pitching velocity
    leadVelLocal=leadVel*rotMat; trailVelLocal=trailVel*rotMat;
    thetDot=(leadVelLocal(3)-trailVelLocal(3))/cordL;
    
    %% get key coeffs
    CL=A*sin(2*(alph-alph0));
    CD=B-C*cos(2*(alph-alph0));
    Crot=pi/2;

    %% get bulk force on this blade
    qCordVelMag=norm(norm(qCordVel));
    FL=.5*rho*spanThick*cordL*CL*qCordVelMag*cross(qCordVel,yDir); % lift
    FD=-.5*rho*spanThick*cordL*CD*qCordVelMag*qCordVel; % drag
    Frot=rho*spanThick*(cordL^2)*thetDot*qCordVelMag*Crot*zDir; % force from pitching velocity
    
    Ftotal=FL+FD+Frot;

    %% distribute forces to elements
    forceFracs=distribFunc(xRel).*As;
    forceFracs=forceFracs/sum(forceFracs);
    elemForces=Ftotal.*forceFracs;
    triangleForces(elemInds,:)=triangleForces(elemInds,:)+elemForces;
end
%% distribute the forces to the masses
unCombinedForces=repmat(triangleForces/3,[3 1]);
indList=[iList;jList;kList];
Forces=splitapply(@(x)sum(x,1),unCombinedForces,indList);

end