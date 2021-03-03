function f=fluxL(u,sigma,model)
% flux left of the discontinuity.
% sigma is a random variable.

% f=u;

oil=lambda(u,sigma);
water=lambda(1-u,sigma);
n=oil+water;

if strcmp(model,'m2') | strcmp(model,'m3')
    kL=1+sigma(2);
    f=(oil./n).*(1-kL.*water);
else
    f=(oil./n).*(1-water);
end
