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
    
    colormap(p.Results.colormap)
    
    dims = size(psi);
    if length(dims) == 4
        for f=1:1:dims(1)
            if p.Results.cutoff == -1
                dens = abs(squeeze(psi(f,:,:,:))).^2;
                dlevel = p.Results.level;
            else
                [~,dens] = calc_psi_cutoff(squeeze(psi(f,:,:,:)),p.Results.cutoff);
                dlevel = p.Results.level*mean(dens(:));
            end
            [faces,v] = isosurface(gridx,gridy,gridz,dens,dlevel);
            q = patch('Faces',faces,'Vertices',v);
            q.EdgeColor = 'none';
            q.FaceAlpha = p.Results.facealpha{f};
            q.FaceColor = p.Results.facecolor{f};
            isonormals(gridx,gridy,gridz,dens,q)
        end
    else
        if p.Results.cutoff == -1
            dens = abs(psi).^2;
            dlevel = p.Results.level;
        else
            [~,dens] = calc_psi_cutoff(psi,p.Results.cutoff);
            dlevel = p.Results.level*mean(dens(:));
        end
        [f,v] = isosurface(gridx,gridy,gridz,dens,dlevel);
        q = patch('Faces',f,'Vertices',v);
        q.EdgeColor = 'none';
        q.FaceAlpha = p.Results.facealpha;
        q.FaceColor = p.Results.facecolor;
        isonormals(gridx,gridy,gridz,dens,q)
    end
    
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
        if length(dims) == 4
            phase = angle(squeeze(psi(p.Results.sideslice,:,:,:)));
        else
            phase = angle(psi(:,:,:));
        end
        h = slice(gridx,gridy,gridz,phase,max(gridx),[],[]);
        set(h,'edgecolor','none');
        set(h,'FaceLighting','none');
        h = slice(gridx,gridy,gridz,phase,[],max(gridy),[]);
        set(h,'edgecolor','none');
        set(h,'FaceLighting','none');
        h = slice(gridx,gridy,gridz,phase,[],[],min(gridz));
        set(h,'edgecolor','none');
        set(h,'FaceLighting','none');
    end
    if (p.Results.sideslice ~= 0 && strcmp(p.Results.slicetype , 'dens' ))
        hold on
        if length(dims) == 4
            Nx = dims(2);
            Ny = dims(3);
            Nz = dims(4);
            dens = abs(squeeze(psi(p.Results.sideslice,:,:,:))).^2;
        else
            Nx = dims(1);
            Ny = dims(2);
            Nz = dims(3);
            dens = abs(psi(:,:,:)).^2;
        end
        h = slice(gridx,gridy,gridz,repmat(dens(:,Ny,:),[1,Nx,1]),max(gridx),[],[]);
        set(h,'edgecolor','none');
        set(h,'FaceLighting','none');
        h = slice(gridx,gridy,gridz,repmat(dens(Nx,:,:),[Ny,1,1]),[],max(gridy),[]);
        set(h,'edgecolor','none');
        set(h,'FaceLighting','none');
        h = slice(gridx,gridy,gridz,repmat(dens(:,:,1),[1,1,Nz]),[],[],min(gridz));
        set(h,'edgecolor','none');
        set(h,'FaceLighting','none');
        caxis([0,max(dens(:))]);
    end
    if (p.Results.sideslice ~= 0 && strcmp(p.Results.slicetype , 'avgdens' ))
        hold on
        if length(dims) == 4
            Nx = dims(2);
            Ny = dims(3);
            Nz = dims(4);
            dens = abs(squeeze(psi(p.Results.sideslice,:,:,:))).^2;
        else
            Nx = dims(1);
            Ny = dims(2);
            Nz = dims(3);
            dens = abs(psi(:,:,:)).^2;
        end
        h = slice(gridx,gridy,gridz,repmat(mean(dens,2)/mean(dens(:)),[1,Nx,1]),max(gridx),[],[]);
        set(h,'edgecolor','none');
        set(h,'FaceLighting','none');
        h = slice(gridx,gridy,gridz,repmat(mean(dens,1)/mean(dens(:)),[Ny,1,1]),[],max(gridy),[]);
        set(h,'edgecolor','none');
        set(h,'FaceLighting','none');
        h = slice(gridx,gridy,gridz,repmat(mean(dens,3)/mean(dens(:)),[1,1,Nz]),[],[],min(gridz));
        set(h,'edgecolor','none');
        set(h,'FaceLighting','none');
        caxis([0,max(dens(:))]);
    end
end