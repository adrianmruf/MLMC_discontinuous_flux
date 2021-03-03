function f=fluxR(u,sigma,model)
% flux right of the discontinuity.
% sigma is a random variable.

% f=0.5*u.^2;

oil=lambda(u,sigma);
water=lambda(1-u,sigma);
n=oil+water;

if strcmp(model,'m2') | strcmp(model,'m3')
    kR=2+sigma(3);
    f=(oil./n).*(1-kR.*water);
else
    f=(oil./n).*(1-2.*water);
end