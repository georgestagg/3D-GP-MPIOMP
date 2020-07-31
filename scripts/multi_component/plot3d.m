i=11;
[~,~,~,psi1,~,~,~,~] = getWF('/data/ngs54/research/lattice_space/1',i,'prefix','imag','magnetic',0,'fnum',1);
[gridx,gridy,gridz,psi2,potential,magx,magy,magz] = getWF('/data/ngs54/research/lattice_space/1',i,'prefix','imag','fnum',2);

plot3d(gridx,gridy,gridz,psi1,'facecolor','red');
hold on
plot3d(gridx,gridy,gridz,psi2,'facecolor','blue');
