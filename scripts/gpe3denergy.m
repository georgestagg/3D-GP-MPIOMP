function [totalE,kinE,potE,interE] = gpe3denergy(gridx,gridy,gridz,psi,pot)
dx=gridx(2)-gridx(1);
dy=gridy(2)-gridy(1);
dz=gridz(2)-gridz(1);

[FX, FY, FZ] = gradient(psi,dx,dy,dz);
dpsi = FX.*conj(FX)+FY.*conj(FY)+FZ.*conj(FZ);
kin=0.5*dpsi;

pot=(pot).*psi.*conj(psi);
inter=0.5*psi.*conj(psi).*psi.*conj(psi);

kinE=trapz(kin(:))*dx*dy*dz;
potE=trapz(pot(:))*dx*dy*dz;
interE=trapz(inter(:))*dx*dy*dz;
totalE=real(kinE+potE+interE);

end 