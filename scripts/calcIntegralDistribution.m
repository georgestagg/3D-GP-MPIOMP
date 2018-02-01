function [Fx,Fk] = calcIntegralDistribution(psi)
Fk=0;
Fx=0;
dims = size(psi);
Nx = dims(1);
Ny = dims(2);
Nz = dims(3);

v1=fftn(psi);
v1=v1/(sqrt(Nx*Ny*Nz));
v1=fftshift(v1);
    
kx=(-Nx/2:(Nx/2-1));
ky=(-Nx/2:(Ny/2-1));
kz=(-Nx/2:(Nz/2-1));

[mkx,mky,mkz] = meshgrid(kx,ky,kz);
k2=mkx.^2+mky.^2+mkz.^2;
for kwv=1:Nx
    w1 = abs(v1(k2 < kwv^2)).^2;
    Fk(kwv+1)=sum(w1(:));
    Fx(kwv+1)=kwv;
end