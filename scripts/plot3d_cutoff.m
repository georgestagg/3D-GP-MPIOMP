function plot3d_cutoff(gridx,gridy,gridz,psi,kc)
dims = size(psi);
Nx = dims(1);
Ny = dims(2);
Nz = dims(3);
v1=fftn(psi);
v1=v1/(sqrt(Nx*Ny*Nz));
v1=ifftshift(v1);
    
kx=(-(Nx-1)/2:(Nx-1)/2);
ky=(-(Ny-1)/2:(Ny-1)/2);
kz=(-(Nz-1)/2:(Nz-1)/2);
for j=1:Nx
    for k=1:Ny
        for l=1:Nz
            k2=kx(j)^2+ky(k)^2+kz(l)^2;
            v1(j,k,l)=v1(j,k,l)*max(1-k2/kc^2,0);
        end
    end  
end
v1=ifftshift(v1);
v1=v1*(sqrt(Nx*Ny*Nz));
u1=ifftn(v1);
ndens=abs(u1).^2;
[f,v] = isosurface(gridx,gridy,gridz,ndens,0.05*mean(ndens(:)));
q = patch('Faces',f,'Vertices',v);
q.EdgeColor = 'none';
q.FaceAlpha = '0.5';
q.FaceColor = 'red';
axis([min(gridy) max(gridy) min(gridx) max(gridx) min(gridz) max(gridz)]);
view([-40 35])
camlight;
lighting gouraud;
daspect([1,1,1]);
end