function u=initialdata(x)
% u is a vector approximating the initial datum by taking cell averages of
% the piecewise constant initial datum u0(x)=0.8 if -0.9<=x<=-0.2,
% u0(x)=0.4 otherwise.

dx=x(2)-x(1);

u=zeros(1,length(x)-1);
for i=1:length(x)-1
    if x(i)<=-0.9 && x(i+1)>-0.9
        u(i)=(0.4*(-0.9-x(i))+0.8*(x(i+1)+0.9))/dx;
    elseif x(i)>-0.9 && x(i+1)<=-0.2
        u(i)=0.8;
    elseif x(i)<=-0.2 && x(i+1)>-0.2
        u(i)=(0.8*(-0.2-x(i))+0.4*(x(i+1)+0.2))/dx;
    else
        u(i)=0.4;
    end
end

