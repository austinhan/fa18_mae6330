fileid=fopen('part.txt');
cell=fscanf(fileid,'%f %f %f %f');
x=cell(1:5:end);
y=cell(2:5:end);
t=cell(5:5:end);
plot(t,y)
