% script to compute an approximation to the mean of the random conservation
% law with discontinuous coefficient using a trapezoidal rule in the
% parameter domain. 
% 'm1' is for the first experiment in the paper with random position of the 
% discontinuity, 'm2' for the second experiment with random absolute
% permeabilities.

model = 'm1';
N=2^12;             % number of grid points in the physical domain
xmin=-1; xmax=1;    % bounary points of physical domain
dx=(xmax-xmin)/N;   % mesh width

x=xmin:dx:xmax;
u0=initialdata(x);  % initialize initial data
T = 0.2;            % final time of compuation
if strcmp(model,'m1')
    Nmc = 200;      % number of grid points in parameter domain            
else
    Nmc = 60;       % number of grid points in parameter domain, in each dimension
end

Ufas = gauss_discclaw(u0,x,T,Nmc,model); % use method to compute approximation of mean by Gaussian quadrature
xfas=x; 
if strcmp(model,'m1')
    save 'testfasitrandposition.mat';     % save data
else
    save 'testfasitrandabsolutepermeability.mat';     % save data
end