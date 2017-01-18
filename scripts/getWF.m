function [gridx,gridy,gridz,psi,potential] = getWF(dirarg,frame)

dirarg = regexprep(dirarg, '/$', '');
datalocation = strcat(dirarg, '/psi.%06d.dat');
fname = sprintf(datalocation,frame);
gridx = ncread(fname,'gx');
gridy = ncread(fname,'gy');
gridz = ncread(fname,'gz');
real = ncread(fname,'real');
imag = ncread(fname,'imag');
potential = ncread(fname,'pot');
psi = real + 1i.*imag;
%permute data to matlab's expected (y,x,z)
psi = permute(psi,[2,1,3]);
potential = permute(potential,[2,1,3]);
fclose('all');
end
