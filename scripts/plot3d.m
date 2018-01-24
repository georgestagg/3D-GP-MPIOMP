function plot3d(gridx,gridy,gridz,psi)
    dens = abs(psi).^2;
    [gx,gy,gz] = meshgrid(gridx,gridy,gridz);
    [f,v] = isosurface(gridx,gridy,gridz,dens,0.5);
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