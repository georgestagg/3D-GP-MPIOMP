function [gridx,gridy,gridz,psi,potential,magx,magy,magz] = getWF(dirarg,frame,varargin)
p = inputParser;
p.KeepUnmatched = true;
addRequired(p,'dirarg');
addRequired(p,'frame');
addParameter(p,'prefix','psi');
addParameter(p,'magnetic','0');
parse(p,dirarg,frame,varargin{:});
dirarg = regexprep(dirarg, '/$', '');
datalocation = strcat(dirarg, '/',p.Results.prefix,'.%06d.nc');
fname = sprintf(datalocation,frame);
gridx = ncread(fname,'gx');
gridy = ncread(fname,'gy');
gridz = ncread(fname,'gz');
real = ncread(fname,'fluid_001_real');
imag = ncread(fname,'fluid_001_imag');
magx=0;
magy=0;
magz=0;
if(p.Results.magnetic ~= 0)
    magx = ncread(fname,'mag_x');
    magy = ncread(fname,'mag_y');
    magz = ncread(fname,'mag_z');
    magx = permute(magx,[2,1,3]);
    magy = permute(magy,[2,1,3]);
    magz = permute(magz,[2,1,3]);
end
potential = ncread(fname,'pot');
psi = real + 1i.*imag;
%permute data to matlab's expected (y,x,z)
psi = permute(psi,[2,1,3]);
psi = squeeze(psi);
potential = permute(potential,[2,1,3]);
fclose('all');
end
