function f=fluxR(u,sigma,model)
% flux right of the discontinuity.
% sigma(1:2) is a random variable.

% f=0.5*u.^2;

oil=lambda(u);
water=lambda(1-u);
n=oil+water;

if strcmp(model,'m2')
    kR=2+sigma(2);
    f=(oil./n).*(1-kR.*water);
else
    f=(oil./n).*(1-2.*water);
end