function saveWF(fname,gridx,gridy,gridz,psi,potential,step,time)
%permute data to fortran's format
psi = permute(psi,[2,1,3]);
potential = permute(potential,[2,1,3]);
dim = size(psi);

nccreate(fname,'step','Datatype','int32');
ncwrite(fname,'step',step);
nccreate(fname,'time');
ncwrite(fname,'time',time);

nccreate(fname,'gx','Dimensions',{'x_dim',dim(1)});
ncwrite(fname,'gx',gridx);
nccreate(fname,'gy','Dimensions',{'y_dim',dim(2)});
ncwrite(fname,'gy',gridy);
nccreate(fname,'gz','Dimensions',{'z_dim',dim(3)});
ncwrite(fname,'gz',gridz);
nccreate(fname,'real','Dimensions',{'x_dim',dim(1),'y_dim',dim(2),'z_dim',dim(3)});
ncwrite(fname,'real',real(psi));
nccreate(fname,'imag','Dimensions',{'x_dim',dim(1),'y_dim',dim(2),'z_dim',dim(3)});
ncwrite(fname,'imag',imag(psi));
nccreate(fname,'pot','Dimensions',{'x_dim',dim(1),'y_dim',dim(2),'z_dim',dim(3)});
ncwrite(fname,'pot',potential);
fclose('all');
end
