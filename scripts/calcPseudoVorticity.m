[grfunction [pv,pvx,pvy,pvz] = calcPseudoVorticity(gridx,psi)
    dx = gridx(5)-gridx(4);
    [FXR,FYR,FZR] = gradient(real(psi),dx);
    [FXI,FYI,FZI] = gradient(imag(psi),dx);
    cp = cross([FXR(:),FYR(:),FZR(:)],[FXI(:),FYI(:),FZI(:)],2);
    pv = reshape(sqrt(sum(cp.^2,2)),size(FXR));
    pvx = reshape(cp(:,1),size(FXR));
    pvy = reshape(cp(:,2),size(FXR));
    pvz = reshape(cp(:,3),size(FXR));
end