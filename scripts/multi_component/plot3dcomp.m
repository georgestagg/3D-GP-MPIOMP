i=33;
[~,~,~,psi1,~,~,~,~] = getWF('/data/ngs54/research/2compsphere/',i,'prefix','imag','magnetic',0,'fnum',1);
[gridx,gridy,gridz,psi2,potential,magx,magy,magz] = getWF('/data/ngs54/research/2compsphere/',i,'prefix','imag','fnum',2);
[mgx,mgy,mgz] = meshgrid(gridx,gridy,gridz);
rx = gridx(end);
psi1(mgx.^2+mgy.^2+mgz.^2 > (rx-0.5).^2) = 1.0;
psi2(mgx.^2+mgy.^2+mgz.^2 > (rx-0.5).^2) = 1.0;
plot3d(gridx,gridy,gridz,psi1,'facecolor','red');
hold on
plot3d(gridx,gridy,gridz,psi2,'facecolor','blue');