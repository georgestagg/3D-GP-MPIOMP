function Lz = calcAngularMomentum(gridx,gridy,gridz,psi)
    dx = gridx(5)-gridx(4);
    [X,Y,~] = meshgrid(gridx,gridy,gridz);
    [FX,FY,~] = gradient(psi,dx);
    Lz = -1i*(conj(psi).*X.*FY - conj(psi).*Y.*FX);
end