function [totalE,kinE,potE,interE] = gpe3denergy_dt(dirarg,startno,stride,endno)
totalE = [];
kinE = [];
potE = [];
for i=startno:stride:endno
    [gridx,gridy,gridz,psi,potential] = getWF(dirarg,i);
    fprintf('read %d\n',i);
    [tE,kE,pE,lE] = gpe3denergy(gridx,gridy,gridz,psi,potential);
    j = (i+(stride-startno))/stride;
    totalE(j) = tE;
    kinE(j) = kE;
    potE(j) = pE;
    interE(j) = lE;
end
fclose('all');
end