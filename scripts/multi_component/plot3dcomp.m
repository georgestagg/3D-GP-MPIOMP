i=1;
[~,~,~,psi1,~] = getWF('/data/ngs54/research/normtest/2',i,'prefix','imag','fnum',1);
[gridx,gridy,gridz,psi2,potential] = getWF('/data/ngs54/research/normtest/2',i,'prefix','imag','fnum',2);
[mgx,mgy,mgz] = meshgrid(gridx,gridy,gridz);
clf
subplot(2,3,1)
plot3d(gridx,gridy,gridz,psi1,'facecolor','red','level',0.02,'facealpha',0.5);
axis image
subplot(2,3,4)
plot3d(gridx,gridy,gridz,psi2,'facecolor','blue','level',0.02,'facealpha',0.5);
axis image
subplot(2,3,2)
imagesc(gridx,gridy,squeeze(angle(psi1(:,:,end/2))))
axis image
subplot(2,3,5)
imagesc(gridx,gridy,squeeze(angle(psi2(:,:,end/2))))
axis image
subplot(2,3,3)
imagesc(gridx,gridy,squeeze(abs(psi1(:,:,end/2))).^2)
axis image
subplot(2,3,6)
imagesc(gridx,gridy,squeeze(abs(psi2(:,:,end/2))).^2)
axis image