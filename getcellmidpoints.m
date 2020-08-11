function xmid=getcellmidpoints(x)
% gives the cell mid points given a discretization of the domain into
% intervals.
% Note that length(xmid)=length(x)-1.

dx=x(2)-x(1);
xmid=x(1)+dx/2:dx:x(end)-dx/2;
