function [u,U]=FVM(u0,x,T,sigma,model)
%  Calculates an approximate solution of u_t + f(k(x),u)_x = 0 with initial
% data u0, up to a time T. x is the space discretization. Periodic boundary
%  conditions are used. sigma is a random variable containing 
% information about the discontinuous coefficient k and the relative
% permeabilities

% 'u' is approximation at final time, 'U' (optional output argument) matrix
% with approximation at all time steps.


N=length(u0);
dx=x(2)-x(1);
% use CFL number 0.2
dt = 0.2*dx;
% number of time steps
Nstep=ceil(T/dt);
dt=T/Nstep;
u=u0;
if nargout>1
    U=zeros(Nstep+1,N);
    U(1,:)=u;
end


% loop over time steps
for i=1:Nstep
    u=explstep(u,x,dt,dx,sigma,model);	% compute approximation at next time step using explicit method
    if nargout>1 % if more than one output argument is given, save all time steps
        U(i+1,:)=u;
    end
end