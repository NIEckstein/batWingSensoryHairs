function [L,pVal,Lsamp] = myLeesL(X,Y,weightFun,numPerms)
% X contains spatial coordinates. 
% size(X) should be [numLocations numSpatialDimensions]
% Y contains independent variable values for each spatial location in X.
% size(Y) should be [numLocations 2]
% weightFun should be a handle to a vectorized function that converts
% distances to spatial connectivity weights
% numPerms specifies how many random spatial permutations to run for
% generating an imperical distribution for the p-value
% L will be Lee's L statistic, pVal will be the p-value 

%% compute weight matrix
distMat=squareform(pdist(X,"euclidean"));
W=weightFun(distMat);
W=W./sum(W,2);
%% compute L for this data
L=Lfun(Y);

%% do spatial permutations to get pvalue
[~,xInds]=sort(rand(size(W,1),1,numPerms),1);
xInds=repmat(xInds,1,2);
yInds=ones(size(xInds)).*[1 2];
linInds=sub2ind(size(Y),xInds,yInds);
Yperms=Y(linInds);
Lsamp=Lfun(Yperms);
pVal=sum(Lsamp(:)>=L)/numPerms;

%% function for computing L
    function LeesL=Lfun(V)
        Vbar=pagemtimes(W,V);
        Vbar=Vbar-mean(Vbar,1);
        LeesL=sum(prod(Vbar,2),1)./sqrt(prod(sum(Vbar.^2,1),2));
    end


end