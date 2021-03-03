% script to compute an approximation to the mean of the random conservation
% law with discontinuous coefficient using the MLMCFVM with 8 levels. 
% 'm1' is for the first experiment in the paper with random position of the 
% discontinuity, 'm2' for the second experiment with random absolute
% permeabilities, and 'm3' for the third experiment where everything is
% random

model = 'm1';
L=8;
coarsestlevel=5;

N=2^(coarsestlevel+L);      % number of grid points in the physical domain
xmin=-1; xmax=1;            % bounary points of physical domain
dx=(xmax-xmin)/N;           % mesh width

xfas=xmin:dx:xmax;

[Ufas,V,V2,x,cputime]=mlmc_method(L,model,coarsestlevel,1);

% save data
if strcmp(model,'m1')
    save 'testfasitrandposition.mat';   
elseif strcmp(model,'m2')
    save 'testfasitrandabsolutepermeability.mat';
else
    save 'testfasitrandeverything.mat';
end
