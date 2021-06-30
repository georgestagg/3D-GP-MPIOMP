function [totalE,kinE,potE,interE,chemPotE] = calcEnergy(gridx,gridy,gridz,psi,varargin)
p = inputParser;
p.KeepUnmatched = true;
addRequired(p,'psi');
addRequired(p,'gridx');
addRequired(p,'gridy');
addRequired(p,'gridz');
addParameter(p,'magnetic',0);
addParameter(p,'potential',0);
addParameter(p,'RHSType',0);
addParameter(p,'N',[0,0,0]);
parse(p,gridx,gridy,gridz,psi,varargin{:});

dx=gridx(2)-gridx(1);
dy=gridy(2)-gridy(1);
dz=gridz(2)-gridz(1);
Lx=max(gridx)-min(gridx);
Ly=max(gridy)-min(gridy);
Lz=max(gridz)-min(gridz);
[mgx,mgy,mgz] = meshgrid(gridx,gridy,gridz);

if(p.Results.RHSType == 0)
    [FX, FY, FZ] = gradient(psi,dx,dy,dz);
    kin = 0.5*FX.*conj(FX)+FY.*conj(FY)+FZ.*conj(FZ);
    pot=p.Results.potential.*psi.*conj(psi);
    inter=0.5*abs(psi.*conj(psi).*psi.*conj(psi));
    chemPot=-(psi.*conj(psi));
end

if(p.Results.RHSType == 1)
    [FX, FY, FZ] = gradient(psi,dx,dy,dz);
    kin = FX.*conj(FX)+FY.*conj(FY)+FZ.*conj(FZ);
    pot=p.Results.potential.*psi.*conj(psi);
    inter=0.5*abs(psi.*conj(psi).*psi.*conj(psi));
    chemPot=-(psi.*conj(psi));
end

if(p.Results.RHSType == 3)
    Ux = exp(-1i.*pi*(p.Results.N(2).*mgz./(Lx*Lz) - p.Results.N(3).*mgy./(Lx*Ly)).*mgx);
    Uy = exp(-1i.*pi*(-p.Results.N(1).*mgz./(Ly*Lz) + p.Results.N(3).*mgx./(Lx*Ly)).*mgy);
    Uz = exp(-1i.*pi*(p.Results.N(1).*mgy./(Ly*Lz) - p.Results.N(2).*mgx./(Lx*Lz)).*mgz);
    [FX, ~, ~] = gradient(Ux.*psi,dx,dy,dz);
    [~, FY, ~] = gradient(Uy.*psi,dx,dy,dz);
    [~, ~, FZ] = gradient(Uz.*psi,dx,dy,dz);
    kin = abs(conj(Ux).*FX).^2+abs(conj(Uy).*FY).^2+abs(conj(Uz).*FZ).^2;
    pot=p.Results.potential.*psi.*conj(psi);
    inter=0.5*(psi.*conj(psi)-1).^2;
    chemPot=0;
end

kinE=trapz(kin(:))*dx*dy*dz;
potE=trapz(pot(:))*dx*dy*dz;
interE=trapz(inter(:))*dx*dy*dz;
chemPotE=trapz(chemPot(:))*dx*dy*dz;
totalE=real(kinE+potE+interE+chemPotE);
end 