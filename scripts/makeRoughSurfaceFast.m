function [K1,K2] = makeRoughSurfaceFast(gridx,gridy,l,sigma)
%% Reference: Kroese, D. P., & Botev, Z. I. (2015). Spatial Process Simulation. DOI: 10.1007/978-3-319-10064-7_12
    function K = Kexpquad(posvec)
        K=sigma^2*exp(-(posvec(1)^2+posvec(2)^2)/(2*l^2));
    end
    m = length(gridx);
    n = length(gridy);
    Rows=zeros(m,n); Cols=Rows;
    for i=1:n
        for j=1:m
            Rows(j,i)=Kexpquad([gridx(i)-gridx(1),gridy(j)-gridy(1)]); 
            Cols(j,i)=Kexpquad([gridx(1)-gridx(i),gridy(j)-gridy(1)]); 
        end
    end

    BlkCirc_row=[Rows, Cols(:,end:-1:2);Cols(end:-1:2,:), Rows(end:-1:2,end:-1:2)];

    lam=real(fft2(BlkCirc_row))/(2*m-1)/(2*n-1);
    if abs(min(lam(lam(:)<0)))>10^-6
        error('Could not find positive definite embedding!')
    else
        lam(lam(:)<0)=0; lam=sqrt(lam);
    end

    F=fft2(lam.*complex(randn(2*m-1,2*n-1),randn(2*m-1,2*n-1)));
    F=F(1:m,1:n);
    K1=real(F); K2=imag(F);
end

