function matplay(x,U,ylimits,T)
% Makes a movie of the matrix U(x,t)
figure(1);
cla;
if nargin<3,
  mx=max(max(U));
  mi=min(min(U));
else
  mx=ylimits(2);
  mi=ylimits(1);
end;
[n,m]=size(U);
nx=length(x);
if (nx==m)
  U=U';
  nplay=n;
else
  nplay=m;
end;
if nargin==4,
    t=linspace(0,T,nplay);
end
axis([x(1) x(nx) mi mx]);
h=gca;
set(h,'Box','on');
lne=line(x,U(:,1));
if nargin==4,
    tit=title(sprintf('T=%2.2f',t(1)));
    set(tit,'Interpreter','Latex','FontSize',20);
end;
for i=2:nplay,
  set(lne,'YData',U(:,i));
  if nargin==4,
      set(tit,'String',sprintf('T=%2.2f',t(i)));
  end
  drawnow;
end;
