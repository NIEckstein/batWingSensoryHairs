function [max,maxLoc,yEst] = fitBezierByLM(x,y,numSteps)
% get an initial guess
h0=mean(y)*10;
absy=abs(y);
p0=sum(absy.*x)/sum(absy);
% [~,idx]=max(abs(y));
% p0=x(idx);
% h0=10*y(idx);
vars=[p0;h0];
lambda=.1*eye(2);
[yEst,J] = bezierFun(vars,x);
oldCost=sum((yEst-y).^2);
for k=1:numSteps
    % % take a regularized gauss newton step and keep it if it worked
    JtransJ=J'*J;
    JtransJ=JtransJ+lambda;% regularize 
    varsNext=vars-JtransJ\(J'*(yEst-y));
    % project to constraint
    if varsNext(1)>.995
        varsNext(1)=.995;
    elseif varsNext(1)<.005
        varsNext(1)=.005;
    end
    [yEstNext,JNext] = bezierFun(varsNext,x);
    newCost=sum((yEstNext-y).^2);
    if newCost<oldCost
        oldCost=newCost;
        vars=varsNext;
        lambda=lambda/3;
        yEst=yEstNext;
        J=JNext;
    else
        lambda=lambda*2;
    end
    
   
end
if nargout>2
yEst = besierFun(vars,x);
end
max=vars(2)/2;
maxLoc=vars(1)/2+.25;
end