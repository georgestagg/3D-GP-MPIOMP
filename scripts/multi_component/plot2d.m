function plot2d(gridx,gridy,psi,varargin)
    
    p = inputParser;
    p.KeepUnmatched = true;
    addRequired(p,'psi');
    addRequired(p,'gridx');
    addRequired(p,'gridy');
    addParameter(p,'z',1);
    addParameter(p,'colormap','parula');
    parse(p,psi,gridx,gridy,varargin{:});
    dims = size(psi);
    figure;
    if length(dims) == 4
        for f=1:1:dims(1)
            subplot(1,dims(1),f);
            dens = abs(squeeze(psi(f,:,:,:))).^2;
            colormap(p.Results.colormap)
            imagesc(gridx,gridy,dens(:,:,p.Results.z));
            axis image
            axis xy
        end
    else
        dens = abs(psi(:,:,:)).^2;
        colormap(p.Results.colormap)
        imagesc(gridx,gridy,dens(:,:,p.Results.z));
        axis image
        axis xy
    end

end