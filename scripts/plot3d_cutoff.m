function ndens = plot3d_cutoff(gridx,gridy,gridz,psi,kc,avg)
    dims = size(psi);
    Nx = dims(1);
    Ny = dims(2);
    Nz = dims(3);
    [~,ndens] = calc_psi_cutoff(psi,kc);
    dens=abs(psi).^2;
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

end