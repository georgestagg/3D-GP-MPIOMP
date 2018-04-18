function grid = addCylinderPotential(grid,gridx,gridy,gridz)
x0 = 32;
y0 = 32;
r=30;
[mgx,mgy,~] = meshgrid(gridx,gridy,gridz);
pot = 50*(sqrt((mgx-x0).^2 + (mgy-y0).^2)>r);
grid(pot>1) = 0;
end