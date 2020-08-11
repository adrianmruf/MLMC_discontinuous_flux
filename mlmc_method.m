function [Umlmc,Var1,Var2,xmlmc,timeused]=mlmc_method(L,model,coarsestlevel,paralleloption)
% Estimate the expectation U and the variance V using the MLMCFVM with L
% levels with 2^coarsestlevel meshpoints at the coarsest level.
% Input parameter model specifies which type of randomess is investigated,
% Implemented are 'm1' and 'm2', 'm1' is random position of the
% discontinuity of the coefficient k and 'm2' is random absolute
% permeabilities to the left and right of the discontinuity located at x=0.
% Input parameter 'paralleloption' is either 0 or 1, default is 0, no
% parallel loop.

if nargin<4
    paralleloption=0;
    if nargin<3
        coarsestlevel = 5;	% #meshpoints at coarsest level is 2^coarsestlevel
        if nargin<2
            model = 'm1';
            fprintf('You did not specify which model you want, I run default random position');
            if nargin<1
                L=5;
                fprintf('You did not specify how many levels, I run L=5. It might take a while >:)');
            end
        end
    end
end

tic;

% initialize ranges for random variables
if strcmp(model,'m2') % random absolute permeabilities
    a=[-0.3 -0.3];
    b=[0.3,0.3];
else if strcmp(model,'m1') % random position of discontinuity
        a=-0.3;
        b=0.3;
    else     % default to model m1
        a=-0.3;
        b=0.3;
        model ='m1';
        fprintf('Your input model parameter has not been implemented yet, I am running the random position model');
    end
end
       
T=0.2;              % final time where solution is computed
finestlevel=coarsestlevel+L;  
N=2^finestlevel;    % number of mesh points at finest level
xmin=-1; xmax=1;     % spatial domain
dxf=(xmax-xmin)/N;

xmlmc=xmin:dxf:xmax;
xmaster=xmlmc;
Umlmc=zeros(size(xmaster)-[0, 1]); % the results are linearly interpolated onto xmlmc
Usquare=Umlmc;  % to estimate the second moment
Var1=Umlmc;     % attempt to estimate variance in a multilevel fashion
dx0=(xmax-xmin)/2^coarsestlevel;    % cell size at coarsest level
samplenumbers=nsamples(dx0,L-1);    % sample numbers for all the levels

samplelevel=finestlevel;

if paralleloption
    for level=1:L-1 % iterate over the levels, start with the finest resolution
        Ml=samplenumbers(level);
        
        % This is only used in case of model m2
        x1master=xmaster(1:2:end);  % grid at next finest level
        %%%%%%%%
        
        Usquaretemp= zeros(size(xmlmc)-[0, 1]);
        utemp = zeros(size(xmlmc)-[0, 1]);
    
        parfor i=1:Ml    % replace by 'for i=1:Ml' if you don't want to run it in parallel
        % (it is possible to restrict the number of workers used by
        % parpool([number]) or in 'Manage cluster profiles' under the menu
        % 'Parallel')
        
            sigma=a+(b-a).*rand(1,2);   % initialize random paramters
            
            if strcmp(model,'m1')
                sigma(2)=sigma(1); % just one random number for random position
                
                dx=xmaster(2)-xmaster(1);
                
                xL=sigma(1):-dx:-1.;
                xL=flip(xL);
                xR=sigma(1)+dx:dx:1.+dx;
                x=[xL,xR];
                
                x1L=sigma(1):-2*dx:-1.;
                x1L=flip(x1L);
                x1R=sigma(1)+2*dx:2*dx:1.+2*dx;
                x1=[x1L, x1R];
            else % model m2
                x=xmaster;
                x1=x1master;
            end
            u0ml=initialdata(x);
            u0ml1=initialdata(x1);
            
            uml=FVM(u0ml,x,T,sigma,model);    % compute an approximation with FVM
            uml1=FVM(u0ml1,x1,T,sigma,model); % same on next finest level
            ul=interp1(getcellmidpoints(x),uml,getcellmidpoints(xmlmc),'linear','extrap');  % interpolate approximation on finest level
            ul1=interp1(getcellmidpoints(x1),uml1,getcellmidpoints(xmlmc),'linear','extrap');
            dUmlmc=(ul-ul1);
            Usquare = Usquare + (ul.^2-ul1.^2)/Ml;  % compute estimator for 2nd moment in a multilevel fashion
            Umlmc=Umlmc+dUmlmc/Ml;        % add new sample to the sum giving an estimator for the mean in the end.
            utemp = utemp+ dUmlmc/Ml; % only samples of this level
            Usquaretemp = Usquaretemp +(ul-ul1).^2/Ml;
        end
        Var1 = Var1+ Usquaretemp-utemp.^2; % compute variance in a multilevel fashion
        clear Usquaretemp;
        clear utemp;
        samplelevel=samplelevel-1;
        xmaster=x1master;
    end

    % compute estimator for coarsest level
    level=L;
    Ml=samplenumbers(level);
    
    Usquaretemp= zeros(size(xmlmc)-[0, 1]);
    utemp = Usquaretemp;

    parfor i=1:Ml % replace by ordinary for loop if you don't want to run things in parallel
        % (it is possible to restrict the number of workers used by
        % parpool([number]) or in 'Manage cluster profiles' under the menu
        % 'Parallel')
        sigma=a+(b-a).*rand(1,2);
        
        if strcmp(model,'m1')
            sigma(2)=sigma(1); % just one random number for random position
            
            dx=xmaster(2)-xmaster(1);
            
            xL=sigma(1):-dx:-1.;
            xL=flip(xL);
            xR=sigma(1)+dx:dx:1.+dx;
            x=[xL,xR];
        else % model m2
            x=xmaster;
        end
        u0ml=initialdata(x);
        
        uml=FVM(u0ml,x,T,sigma,model);
        ul=interp1(getcellmidpoints(x),uml,getcellmidpoints(xmlmc),'linear','extrap');
        dUmlmc=ul;
        Umlmc=Umlmc+dUmlmc/Ml;    % add new sample normalized
        Usquare = Usquare + ul.^2/Ml;   
        utemp = utemp+ dUmlmc/Ml; % only samples of this level
        Usquaretemp = Usquaretemp +(ul).^2/Ml;
    end
else
    displaylevel=progress(' Level: ');
    for level=1:L-1 % iterate over the levels, start with the finest resolution
        Ml=samplenumbers(level);
        
        x1master=xmaster(1:2:end);
        showprogress(displaylevel,L-level);
        computinglevel=progress(' Computing: ');
        
        Usquaretemp= zeros(size(xmlmc)-[0,1]);
        utemp = zeros(size(xmlmc)-[0,1]);
    
        for i=1:Ml        
            sigma=a+(b-a).*rand(1,2);   % initialize random paramters
            if strcmp(model,'m1')
                sigma(2)=sigma(1); % just one random number for random position
                
                dx=xmaster(2)-xmaster(1);
                
                xL=sigma(1):-dx:-1.;
                xL=flip(xL);
                xR=sigma(1)+dx:dx:1.+dx;
                x=[xL,xR];
                
                x1L=sigma(1):-2*dx:-1.;
                x1L=flip(x1L);
                x1R=sigma(1)+2*dx:2*dx:1.+2*dx;
                x1=[x1L, x1R];
            else % model m2
                x=xmaster;
                x1=x1master;
            end
            u0ml=initialdata(x);
            u0ml1=initialdata(x1);
            showprogress(computinglevel,ceil(100*i/Ml));
            
            uml=FVM(u0ml,x,T,sigma,model);    % compute an approximation with FVM
            uml1=FVM(u0ml1,x1,T,sigma,model); % same on next finest level
            ul=interp1(getcellmidpoints(x),uml,getcellmidpoints(xmlmc),'linear','extrap');  % interpolate approximation on finest level
            ul1=interp1(getcellmidpoints(x1),uml1,getcellmidpoints(xmlmc),'linear','extrap');
            dUmlmc=(ul-ul1);
            Usquare = Usquare + (ul.^2-ul1.^2)/Ml;  % compute estimator for 2nd moment in a multilevel fashion
            Umlmc=Umlmc+dUmlmc/Ml;        % add new sample to the sum giving an estimator for the mean in the end.
            utemp = utemp+ dUmlmc/Ml; % only samples of this level
            Usquaretemp = Usquaretemp +(ul-ul1).^2/Ml;
        end
        Var1 = Var1+ Usquaretemp-utemp.^2; % compute variance in a multilevel fashion
        clear Usquaretemp;
        clear utemp;
        deleteprogress(computinglevel);
        samplelevel=samplelevel-1;
        xmaster=x1master;
    
    end

    % compute estimator for coarsest level
    level=L;
    Ml=samplenumbers(level);
    showprogress(displaylevel,0);
    computinglevel=progress(' Computing: ');
    Usquaretemp= zeros(size(xmlmc)-[0,1]);
    utemp = Usquaretemp;

    for i=1:Ml 
        sigma=a+(b-a).*rand(1,2);
        if strcmp(model,'m1')
            sigma(2)=sigma(1); % just one random number for random position
            
            dx=xmaster(2)-xmaster(1);
            
            xL=sigma(1):-dx:-1.;
            xL=flip(xL);
            xR=sigma(1)+dx:dx:1.+dx;
            x=[xL,xR];
        else % model m2
            x=xmaster;
        end
        u0ml=initialdata(x);
        
        showprogress(computinglevel,ceil(100*i/Ml));
        uml=FVM(u0ml,x,T,sigma,model);
        ul=interp1(getcellmidpoints(x),uml,getcellmidpoints(xmlmc),'linear','extrap');
        dUmlmc=ul;
        Umlmc=Umlmc+dUmlmc/Ml;    % add new sample normalized
        Usquare = Usquare + ul.^2/Ml;
          
        utemp = utemp+ dUmlmc/Ml; % only samples of this level
        Usquaretemp = Usquaretemp +(ul).^2/Ml;
    end
    deleteprogress(computinglevel);
    deleteprogress(displaylevel);

end
    
xmlmc=getcellmidpoints(xmlmc);

Var1 = Var1+ Usquaretemp-utemp.^2;
Var2 = Usquare -Umlmc.^2;
timeused=toc;


% .............. Routines for the display of progress ...............

function d=progress(s) % to print s 0
d.string=s;
fprintf(d.string);
d.number=sprintf(' %3d',0);
fprintf(d.number);


function showprogress(d,x) % to update s x

slett(d.number);
d.number=sprintf(' %3d',x);
fprintf(d.number);


function deleteprogress(d) % to delete s x
slett(d.number);
slett(d.string);

function slett(s)
for i=1:length(s),
    fprintf('\b');
end;

