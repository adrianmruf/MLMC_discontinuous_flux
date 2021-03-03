function u=explstep(u,x,dt,dx,sigma,model)
%
% One step explicit scheme for scheme v=u - dt* Dm f_j
% periodic boundary. Upwing method, assumes f'>0.
%
% sigma - vector of random numbers

N=length(u);

lambda=dt/dx;

J=1;
if strcmp(model,'m2')
    % J=int32(N/2);
    J=find(x==0.);
else
    J=find(x==sigma(1));
end

g=fluxL(u,sigma,model);
f=fluxR(u,sigma,model);

% periodic boundary conditions
IL=1:J-1;
ILm=[N,1:J-2];

uL=u(IL)-lambda*(g(IL)-g(ILm));

uM=fluxRinverse(fluxL(uL(end),sigma,model),sigma,model);

IR=J+1:N;
IRm=J:N-1;
uR=u(IR)-lambda*(f(IR)-f(IRm));

u=[uL,uM,uR];
