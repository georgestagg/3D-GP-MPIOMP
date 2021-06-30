function [totalE,kinE,potE,interE,chemPotE] = calcEnergy_dt(dirarg,startno,stride,endno,varargin)
totalE = [];
kinE = [];
potE = [];
for i=startno:stride:endno
    [gridx,gridy,gridz,psi,~] = getWF(dirarg,i,varargin{:});
    fprintf('read %d\n',i);
    [tE,kE,pE,iE,cE] = calcEnergy(gridx,gridy,gridz,psi,varargin{:});
    j = (i+(stride-startno))/stride;
    totalE(j) = tE;
    kinE(j) = kE;
    potE(j) = pE;
    interE(j) = iE;
    chemPotE(j) = cE;
end
fclose('all');
end