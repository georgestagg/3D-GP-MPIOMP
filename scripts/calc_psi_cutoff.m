function [npsi,ndens] = calc_psi_cutoff(psi,kc)
    dims = size(psi);
    Nx = dims(1);
    Ny = dims(2);
    Nz = dims(3);
    kx=(-Nx/2:(Nx/2-1));
    ky=(-Ny/2:(Ny/2-1));
    kz=(-Nz/2:(Nz/2-1));
    v1=fftshift(fftn(psi));
    [mkx,mky,mkz] = meshgrid(kx,ky,kz);
    k2=mkx.^2+mky.^2+mkz.^2;
    v1=v1.*max(1-k2./kc.^2,0);
    npsi=ifftn(ifftshift(v1));
    ndens=abs(npsi).^2;
end