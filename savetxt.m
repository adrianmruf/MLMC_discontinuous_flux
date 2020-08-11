function savetxt(x,u,filename)
% writes the solution (x,u) as a csv file

fid = fopen(filename,'wt');

for i=1:length(x)
    fprintf(fid,'%.8f %.8f \n',[x(i),u(i)]);
end
fclose(fid);