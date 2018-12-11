fileid=fopen('MASS.TXT');
A=fscanf(fileid,'%f');
fclose(fileid);
plot(A)