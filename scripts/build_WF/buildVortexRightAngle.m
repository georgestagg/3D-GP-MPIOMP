function gridRA = buildVortexRightAngle(gridx,gridy,gridz)

br=64;
z0=0;
x0=-32;
y0=32;
rn=60;

[mgx,mgy,mgz] = meshgrid(gridx,gridy,gridz);
s = ((mgx-x0).^rn+(mgz-z0).^rn).^(1/rn);
d1 = sqrt((mgy-y0).^2+(s+br).^2);
d2 = sqrt((mgy-y0).^2+(s-br).^2);
rr1=get_rr(d1);
rr2=get_rr(d2);
gridRA = rr1.*(mgy-y0+1i.*(s+br)).*rr2.*(mgy-y0+1i.*(s-br));

    function rr = get_rr(r)
        rr=sqrt((0.3437+0.0286.*r.^2)./(1+(0.3333.*r.^2)+(0.0286.*r.^4)));
    end
end