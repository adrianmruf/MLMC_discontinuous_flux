function umean=gauss_discclaw(u0,x,T,Ngrid,model)
% function to compute the mean by a Gaussian quadrature in the parameter
% domain
% replace 'parfor' loops by ordinary 'for' loops if you don't want to run
% the program in parallel.
% model 'm1' has one random parameter, model 'm2' two, hence it computes 
% approximations at Ngrid^2 points.

tic
umean=zeros(size(u0));
if strcmp(model,'m1')
    a = [-0.3,0.3];                     % min and max position of discontinuity
    a = linspace(a(1),a(end),Ngrid);    % vector with Gauss points for positions 
    parfor i=1:Ngrid,
        sigma=[a(i) a(i)]; % position of discontinuity
        
        dx=x(2)-x(1);
        xsigL=a(i):-dx:-1.;
        xsigL=flip(xsigL);
        xsigR=a(i)+dx:dx:1.;
        xsig=[xsigL, xsigR];
        
        u0sig=initialdata(xsig)
        
        [u,~]=FVM(u0sig,xsig,T,sigma,model);% compute approximation to pde for this position
        u=interp1(getcellmidpoints(xsig),u,getcellmidpoints(x),'linear','extrap');
        
        w=weights1d(i,Ngrid);           % get weight of this Gauss point
        umean=umean+w*u;                % add point to sum which approximates E[u]
    end
    wsum=2*Ngrid-2;     
    umean=umean/wsum;                   % to normalize weights correctly
  
else if strcmp(model,'m2')
        a=[-0.3 0.3]; b=[-0.3 0.3];     % min and max absolute permeabilities
        a=linspace(a(1),a(end),Ngrid);  % discretize interval
        b=linspace(b(1),b(end),Ngrid);
        wsum=0;
        for i=1:Ngrid
            parfor j=1:Ngrid
                sigma=[a(i) b(j)];      
                [u,~]=FVM(u0,x,T,sigma,model); % compute approximation
                w=weights(i,j,Ngrid);      % get Gauss weights 
                wsum=wsum+w;               
                umean=umean+w*u;           % add point to sum which gives the mean  
            end
        end
        umean=umean/wsum;                  % normalize weights 

    else
        fprintf('Specify which model!')
    end   
end
toc
