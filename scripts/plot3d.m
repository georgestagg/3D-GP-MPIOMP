function plot3d(gridx,gridy,gridz,psi,lvl)
    dens = abs(psi).^2;
    [f,v] = isosurface(gridx,gridy,gridz,dens,lvl);
    q = patch('Faces',f,'Vertices',v);
    q.EdgeColor = 'none';
    q.FaceAlpha = '1.0';
    q.FaceColor = 'red';
    axis([min(gridy) max(gridy) min(gridx) max(gridx) min(gridz) max(gridz)]);
    view([-40 35])
    camlight;
    lighting gouraud;
    daspect([1,1,1]);
end