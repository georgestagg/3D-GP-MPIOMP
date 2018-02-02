function makeMovie(dirarg,startno,stride,endno)
    dirarg = regexprep(dirarg, '/$', ''); 
    pngfolder = strcat(dirarg, '/png');
    mkdir(pngfolder);
    for i=startno:stride:endno
        try
            [gridx,gridy,gridz,psi,~] = getWF(dirarg,i);
        catch e
            disp(e);
            continue;
        end
        fprintf('read %d\n',i);
        j = i/stride;
        h=figure('visible','off');
        %plotStarandPotpins(gridx,gridy,gridz,psi,potential);
        plot3d_cutoff(gridx,gridy,gridz,psi,4,0.05);
        %plot3d(gridx,gridy,gridz,psi,0.25);
        filename = strcat(pngfolder, '/d%04d.png');
        finalfname = sprintf(filename,j);
        print (h,'-dpng','-r300',finalfname);
        close(h);
    end
end
