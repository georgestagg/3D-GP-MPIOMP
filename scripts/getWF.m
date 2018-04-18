function [gridx,gridy,gridz,psi,potential] = getWF(dirarg,frame,varargin)
p = inputParser;
p.KeepUnmatched = true;
addRequired(p,'dirarg');
addRequired(p,'frame');
addParameter(p,'prefix','psi');
parse(p,dirarg,frame,varargin{:});
dirarg = regexprep(dirarg, '/$', '');
datalocation = strcat(dirarg, '/',p.Results.prefix,'.%06d.nc');
fname = sprintf(datalocation,frame);
gridx = ncread(fname,'gx');
gridy = ncread(fname,'gy');
gridz = ncread(fname,'gz');
real = ncread(fname,'real');
imag = ncread(fname,'imag');
potential = ncread(fname,'pot');
psi = real + 1i.*imag;
%permute data to matlab's expected (y,x,z)
psi = permute(psi,[2,1,3]);
potential = permute(potential,[2,1,3]);
fclose('all');
end
