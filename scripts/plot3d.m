function plot3d(gridx,gridy,gridz,psi,varargin)
    
    p = inputParser;
    p.KeepUnmatched = true;
    addRequired(p,'psi');
    addRequired(p,'gridx');
    addRequired(p,'gridy');
    addRequired(p,'gridz');
    addParameter(p,'level',0.25);
    addParameter(p,'facecolor','red');
    addParameter(p,'facealpha',1);
    addParameter(p,'colormap','parula');
    addParameter(p,'colordata',0);
    addParameter(p,'view',[-55 18]);
    addParameter(p,'sideslice',0);
    addParameter(p,'slicetype','phase');
    addParameter(p,'cutoff',-1);
    parse(p,psi,gridx,gridy,gridz,varargin{:});
    
    if p.Results.cutoff == -1
        dens = abs(psi).^2;
        dlevel = p.Results.level;
    else
        [~,dens] = calc_psi_cutoff(psi,p.Results.cutoff);
        dlevel = p.Results.level*mean(dens(:));
    end
    
    colormap(p.Results.colormap)
    [f,v] = isosurface(gridx,gridy,gridz,dens,dlevel);
    q = patch('Faces',f,'Vertices',v);
    q.EdgeColor = 'none';
    q.FaceAlpha = p.Results.facealpha;
    q.FaceColor = p.Results.facecolor;
    isonormals(gridx,gridy,gridz,dens,q)
    
    if(p.Results.colordata ~= 0)
        isocolors(gridx,gridy,gridz,p.Results.colordata,q)
    end
    
    axis([min(gridx) max(gridx) min(gridy) max(gridy) min(gridz) max(gridz)]);
    daspect([1,1,1]);
    view(p.Results.view)
    camlight;
    lighting gouraud;
    
    if (p.Results.sideslice ~= 0 && strcmp(p.Results.slicetype , 'phase' ))
        hold on
        h = slice(gridx,gridy,gridz,angle(psi),max(gridx),[],[]);
        set(h,'edgecolor','none');
        set(h,'FaceLighting','none');
        h = slice(gridx,gridy,gridz,angle(psi),[],max(gridy),[]);
        set(h,'edgecolor','none');
        set(h,'FaceLighting','none');
        h = slice(gridx,gridy,gridz,angle(psi),[],[],min(gridz));
        set(h,'edgecolor','none');
        set(h,'FaceLighting','none');
    end
    if (p.Results.sideslice ~= 0 && strcmp(p.Results.slicetype , 'dens' ))
        hold on
        dims = size(psi);
        Nx = dims(1);
        Ny = dims(2);
        Nz = dims(3);
        dens = abs(psi).^2;
        h = slice(gridx,gridy,gridz,repmat(dens(:,Ny,:),[1,Nx,1]),max(gridx),[],[]);
        set(h,'edgecolor','none');
        set(h,'FaceLighting','none');
        h = slice(gridx,gridy,gridz,repmat(dens(Nx,:,:),[Ny,1,1]),[],max(gridy),[]);
        set(h,'edgecolor','none');
        set(h,'FaceLighting','none');
        h = slice(gridx,gridy,gridz,repmat(dens(:,:,1),[1,1,Nz]),[],[],min(gridz));
        set(h,'edgecolor','none');
        set(h,'FaceLighting','none');
    if (p.Results.sideslice ~= 0 && strcmp(p.Results.slicetype , 'avgdens' ))
        hold on
        dims = size(psi);
        Nx = dims(1);
        Ny = dims(2);
        Nz = dims(3);
        dens = abs(psi).^2;
        h = slice(gridx,gridy,gridz,repmat(mean(dens,2)/mean(dens(:)),[1,Nx,1]),max(gridx),[],[]);
        set(h,'edgecolor','none');
        set(h,'FaceLighting','none');
        h = slice(gridx,gridy,gridz,repmat(mean(dens,1)/mean(dens(:)),[Ny,1,1]),[],max(gridy),[]);
        set(h,'edgecolor','none');
        set(h,'FaceLighting','none');
        h = slice(gridx,gridy,gridz,repmat(mean(dens,3)/mean(dens(:)),[1,1,Nz]),[],[],min(gridz));
        set(h,'edgecolor','none');
        set(h,'FaceLighting','none');
    end
end