function plotStarandPotpins(gridx,gridy,gridz,psi,potential)
    dens = abs(psi).^2;
    [gx,gy,gz] = meshgrid(gridx,gridy,gridz);
    holes = potential - ((gx-12.75).^2+(gy-12.75).^2+(gz-12.75).^2)/2;
    holes(((gx-12.75).^2+(gy-12.75).^2+(1.7*(gz-12.75)).^2) > 56) = 1;
    p = patch(isosurface(gridx,gridy,gridz,holes,50));
    p.EdgeColor = 'none';
    p.FaceAlpha = '1.0';
    p.FaceColor = 'yellow';
    [f,v] = isosurface(gridx,gridy,gridz,dens,0.002);
    q = patch('Faces',f,'Vertices',v);
    q.EdgeColor = 'none';
    q.FaceAlpha = '0.5';
    q.FaceColor = 'red';
    view([-40 35])
    camlight;
    daspect([1,1,1]);
    axis([4 21 4 21 4 21])
end