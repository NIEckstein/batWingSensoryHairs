function [yEst,errorJacobian] = bezierFun(vars,x)

p=vars(1); h=vars(2);

temp1=(1-2*p)*x;
temp2=sqrt(p^2 + temp1);
t=(-p+temp2)/(1-2*p);
yEst=2*h*t.*(1-t);
%error=yEst-y;

if nargout>1
    dtdp= (p+temp1-temp2)./temp2/(2*p-1)^2;
    dyEstdp=2*h*(1-2*t).*dtdp;
    dyEstdh=yEst/h;
    errorJacobian=[dyEstdp dyEstdh];
end


end