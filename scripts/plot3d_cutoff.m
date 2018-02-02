function ndens = plot3d_cutoff(gridx,gridy,gridz,psi,kc,avg)
dims = size(psi);
Nx = dims(1);
Ny = dims(2);
Nz = dims(3);
dx = gridx(2)-gridx(1);
kx=(-Nx/2:(Nx/2-1));
ky=(-Ny/2:(Ny/2-1));
kz=(-Nz/2:(Nz/2-1));
dens=abs(psi).^2;

v1=fftshift(fftn(psi));
[mkx,mky,mkz] = meshgrid(kx,ky,kz);
k2=mkx.^2+mky.^2+mkz.^2;
v1=v1.*max(1-k2./kc.^2,0);
ndens=abs(ifftn(ifftshift(v1))).^2;

[f,v] = isosurface(gridx,gridy,gridz,ndens,avg*mean(ndens(:)));
q = patch('Faces',f,'Vertices',v);
isonormals(gridx,gridy,gridz,ndens,q);
q.EdgeColor = 'none';
q.FaceAlpha = '1.0';
q.FaceColor = 'red';

view([-55 18])
camlight;
lighting gouraud;

hold on
h = slice(gridx,gridy,gridz,repmat(mean(dens,2)/mean(dens(:)),[1,Nx,1]),max(gridx),[],[]);
set(h,'edgecolor','none');
set(h,'FaceLighting','none');
h = slice(gridx,gridy,gridz,repmat(mean(dens,1)/mean(dens(:)),[Ny,1,1]),[],max(gridy),[]);
set(h,'edgecolor','none');
set(h,'FaceLighting','none');
h = slice(gridx,gridy,gridz,repmat(mean(dens,3)/mean(dens(:)),[1,1,Nz]),[],[],min(gridz));
set(h,'edgecolor','none');
set(h,'FaceLighting','none');

axis([min(gridy) max(gridy) min(gridx) max(gridx) min(gridz) max(gridz)]);
daspect([1,1,1]);
%caxis([0.5,1.7]);

end