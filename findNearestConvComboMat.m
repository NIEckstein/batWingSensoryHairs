function [wOpt] = findNearestConvComboMat(w)
[rDim, cDim]=size(w); numEl=rDim*cDim;
costFun=@(wVar) sum((wVar-w(:)).^2);
LB=zeros(numEl,1); UB=LB+1;
Aeq=repmat(eye(rDim),[1 cDim]); Beq=ones(rDim,1);

wOpt=fmincon(costFun,w(:),[],[],Aeq,Beq,LB,UB);
wOpt = reshape(wOpt, rDim, cDim);



end