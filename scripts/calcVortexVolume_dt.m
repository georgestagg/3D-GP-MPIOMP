function vvdt = calcVortexVolume_dt(dirarg,startno,stride,endno)
vvdt = [];
for i=startno:stride:endno
    [gridx,gridy,gridz,psi,~] = getWF(dirarg,i);
    fprintf('read %d\n',i);
    vv = calcVortexVolume_cutoff(gridx,gridy,gridz,psi,7,0.05);
    j = (i+(stride-startno))/stride;
    vvdt(j) = vv;
end
fclose('all');
end