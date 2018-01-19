function plot3d(gridx,gridy,gridz,psi)
    dens = abs(psi).^2;
    avg = mean(dens(:));
    [gx,gy,gz] = meshgrid(gridx,gridy,gridz);
    [f,v] = isosurface(gridx,gridy,gridz,dens,0.05*avg);
    q = patch('Faces',f,'Vertices',v);
    q.EdgeColor = 'none';
    q.FaceAlpha = '0.5';
    q.FaceColor = 'red';
    view([-40 35])
    camlight;
    daspect([1,1,1]);
end