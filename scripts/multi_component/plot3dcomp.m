i=52;
[~,~,~,psi1,~,~,~,~] = getWF('/data/ngs54/research/2compsphere/',i,'prefix','imag','magnetic',0,'fnum',1);
[gridx,gridy,gridz,psi2,potential,magx,magy,magz] = getWF('/data/ngs54/research/2compsphere/',i,'prefix','imag','fnum',2);
[mgx,mgy,mgz] = meshgrid(gridx,gridy,gridz);
tx = gridx(end)/2.0;
psi1((mgx-tx).^2+(mgy-tx).^2+(mgz-tx).^2 > (tx-0.5).^2) = 1.0;
psi2((mgx-tx).^2+(mgy-tx).^2+(mgz-tx).^2 > (tx-0.5).^2) = 1.0;
plot3d(gridx,gridy,gridz,psi1,'facecolor','red');
hold on
plot3d(gridx,gridy,gridz,psi2,'facecolor','blue');
