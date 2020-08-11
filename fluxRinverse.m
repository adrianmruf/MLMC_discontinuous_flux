function u=fluxRinverse(y,sigma,model)
% approximation of the inverse of the flux to the right of the discontinuity

fun=@(x) fluxR(x,sigma,model) -y;
interval=[0.,1.];
u=fzero(fun,interval);