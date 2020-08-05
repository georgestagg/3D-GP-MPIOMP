i=1;
[~,~,~,psi1] = getWF('/data/ngs54/research/2compsphere',i,'prefix','imag','fnum',1);
%[gridx,gridy,gridz,psi2,potential] = getWF('/data/ngs54/research/2compsphere',i,'prefix','imag','fnum',2);
[mgx,mgy,mgz] = meshgrid(gridx,gridy,gridz);
%rx = gridx(end);
%psi1(mgx.^2+mgy.^2+mgz.^2 > (rx-0.2).^2) = 1.0;
%psi2(mgx.^2+mgy.^2+mgz.^2 > (rx-0.2).^2) = 1.0;
clf
subplot(1,3,1)
plot3d(gridx,gridy,gridz,psi1,'facecolor','red','level',0.04);
axis image
hold on
plot3d(gridx,gridy,gridz,psi2,'facecolor','blue','level',0.04);
axis image
subplot(2,3,2)
imagesc(squeeze(angle(psi1(:,:,1))))
axis image
subplot(2,3,5)
imagesc(squeeze(angle(psi2(:,:,1))))
axis image
subplot(2,3,3)
imagesc(squeeze(abs(psi1(:,:,1))).^2)
axis image
subplot(2,3,6)
imagesc(squeeze(abs(psi2(:,:,1))).^2)
axis image