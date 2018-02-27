function vv = calcVortexVolume_cutoff(gridx,gridy,gridz,psi,kc,avg)
    dx = gridx(2)-gridx(1);
    dims = size(psi);
    Nx = dims(1);
    Ny = dims(2);
    Nz = dims(3);
    [~,ndens] = calc_psi_cutoff(psi,kc);
    vv = trapz(ndens(:)<avg*mean(ndens(:)))*dx*dx*dx;
end