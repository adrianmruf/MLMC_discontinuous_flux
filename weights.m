function w=weights(i,j,n)

% computes (unnormalized) weights for Gaussian quadrature with trapezoidal rule in each grid cell, 2d

w=(i>1)+(i<n)+(j>1)+(j<n);
if w==2 
    w=1;
end
if w==3
    w=2;
end
