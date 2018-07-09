function saveSURF(fname,gridx,gridy,surf)
%permute data to fortran's format
surf = surf';
dim = size(surf);

nccreate(fname,'gx','Dimensions',{'x_dim',dim(1)});
ncwrite(fname,'gx',gridx);
nccreate(fname,'gy','Dimensions',{'y_dim',dim(2)});
ncwrite(fname,'gy',gridy);
nccreate(fname,'surf','Dimensions',{'x_dim',dim(1),'y_dim',dim(2)});
ncwrite(fname,'surf',surf);
fclose('all');
end
