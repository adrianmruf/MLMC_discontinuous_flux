function f=fluxL(u,sigma,model)
% flux left of the discontinuity.
% sigma(1:2) is a random variable.

% f=u;

oil=lambda(u);
water=lambda(1-u);
n=oil+water;

if strcmp(model,'m2')
    kL=1+sigma(1);
    f=(oil./n).*(1-kL.*water);
else
    f=(oil./n).*(1-water);
end
