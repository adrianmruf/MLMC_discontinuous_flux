function w=weights1d(i,n)

% computes weights for Gaussian quadrature with trapezoidal rule in each grid cell, 1d, unnormalized weights

w=(i>1)+(i<n);
% if w==2 
%     w=1;
% end
% if w==3
%     w=2;
% end
