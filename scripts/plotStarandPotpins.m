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
    xv = v(:,1);
    yv = v(:,2);
    zv = v(:,3);
    xv(sqrt((xv-12.75).^2+(yv-12.75).^2+(1.7*(zv-12.75)).^2) > 7) = nan;
    yv(sqrt((xv-12.75).^2+(yv-12.75).^2+(1.7*(zv-12.75)).^2) > 7) = nan;
    zv(sqrt((xv-12.75).^2+(yv-12.75).^2+(1.7*(zv-12.75)).^2) > 7) = nan;
    v(:,1) = xv;
    v(:,2) = yv;
    v(:,3) = zv;
    q = patch('Faces',f,'Vertices',v);
    q.EdgeColor = 'none';
    q.FaceAlpha = '0.75';
    q.FaceColor = 'red';
    view([-40 35])
    camlight;
    daspect([1,1,1]);
    axis([4 21 4 21 4 21])
end