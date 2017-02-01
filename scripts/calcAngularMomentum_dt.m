function angm_dt = calcAngularMomentum_dt(dirarg,startno,stride,endno)
    dirarg = regexprep(dirarg, '/$', ''); 
    for i=startno:stride:endno
        [gridx,gridy,gridz,psi,~] = getWF(dirarg,i);
        dx = gridx(5)-gridx(4);
        fprintf('read %d\n',i);
        j = (i - startno)/stride  + 1;
        angm = calcAngularMomentum(gridx,gridy,gridz,psi);
        angm_dt(j) = trapz(angm(:))*dx*dx*dx;
    end
end
