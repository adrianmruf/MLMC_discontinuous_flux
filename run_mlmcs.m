% Calculates the average of ntry mlmc solutions at finest level i, 
% for i=1,...,maxtest
% 'm1' is the model with the random position of the discontinuity of k, 
% 'm2' is for the model with random absolute permeabilities.

model='m1';             % pick a model, 'm1', 'm2' are implemented
paralleloption = 1;     % if set =1, the inner loop of MLMC_method is run in parallel, if set =0, no parallelization is used.
maxtest=6;  % we test up to maxtest levels. If more than ca. 6 are used one has to 
            % produce a reference solution with higher resolution

ntry=30;    % number of MLMC estimators that are computed at each level. We take
            % the average of the relative errors at each level.
coarsestmesh = 5; % mesh resolution at the coarsest level is 2^coarsestmesh            

rng(1234);  % sets the random generator to use a specific sequence

Rk=zeros(1,maxtest);    % initializing variables for relative error, variance, cpu time
vRk=zeros(1,maxtest);
timeused=zeros(1,maxtest);

% load reference solutions
if strcmp(model,'m1')
    load('testfasitrandposition.mat');  % make sure folder hierarchy is correct
else
    load('testfasitrandabsolutepermeability.mat');
end

for i=1:maxtest 
    R=0;                % temporary variables for relative error and variance of errors.
    var=0;
    for k=1:ntry,
        [U,V,V2,x,tu]=mlmc_method(i,model,coarsestmesh,paralleloption);   % compute an mlmc estimator with i levels
        timeused(i)=timeused(i)+tu;
        Uf=interp1(x,U,getcellmidpoints(xfas),'linear','extrap');
        err=Uf-Ufas;
        % estimating mean by average instead of RMS
        h=R;
        e=100*norm(err,1)./norm(Ufas,1);
        R=((k-1)*R+e)/k;
        var=((k-1)/k)*var+((k-1)/k^2)*(e-h)^2+(1/k)*(e-R)^2;
    end
    timeused(i)=timeused(i)/ntry;
    Rk(i)=R;
    vRk(i)=var;
    %figure; plot(x,U,xfas,Ufas)
    fprintf('\n L=%d finished \n',i);
end
fprintf('\n');
resolutions = 2.^[-coarsestmesh:-1:-(coarsestmesh+maxtest-1)];
rate1 = polyfit(log2(resolutions),log2(Rk),1); % use polyfit to compute the rate error vs. the mesh resolution
ratemesh = abs(rate1(1));
rate2 = polyfit(log2(timeused),log2(Rk),1);  % use polyfit to compute the rate error vs. cputime used
ratework = abs(rate2(1));
if strcmp(model,'m1')
    save('RelativeErrorsfilerandposition','Rk','vRk','timeused','ratemesh', 'ratework');   % save errors, variance and cpu time into a .mat file.
else
    save('RelativeErrorsfilerandabsolutepermeability','Rk','vRk','timeused','ratemesh', 'ratework');   % save errors, variance and cpu time into a .mat file.
end