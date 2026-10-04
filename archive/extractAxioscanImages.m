function  extractAxioscanImages(basepath, varargin)

datapaths = dir(basepath);
datapaths(~[datapaths(:).isdir])=[];
resmax = 2^8-1;
for pathNum = 4:-1:3
    datapath = datapaths(pathNum);
    fprintf('Path number: %d, folder name %s.\n',pathNum,datapath.name);
    if datapath.isdir
        savefolder    = fullfile(datapath.folder,datapath.name,'extracted_cy3');
        datapath      = fullfile(datapath.folder,datapath.name);
    else
        savefolder    = fullfile(datapath.folder,'extracted_cy3');
        datapath      = datapath.folder;
    end
    %-------------------------------------------------------------------------
    imagefiles  = dir(fullfile(datapath, '**', '*.czi'));
    
    if ~isempty(imagefiles)
        for ifile = 1:numel(imagefiles)
            dpfile  = fullfile(imagefiles(ifile).folder, imagefiles(ifile).name);
            [~, nameuse, ~] = fileparts(dpfile);

            dataim  = BioformatsImage(dpfile);
            %select relevant series
            Nseries = dataim.seriesCount;
            allwh   = nan(Nseries, 2);
            allpx   = nan(Nseries, 2);
            for is = 1:Nseries
                dataim.series = is;
                allwh(is, :) = [dataim.width dataim.height];
                allpx(is, :) = dataim.pxSize;
            end
            irel = [1; find(diff(allwh(:,1))>0)+1];
            medw = median(allwh(irel, 1));
            irem = abs((allwh(irel,1) - medw)/medw)>0.3;
            irel(irem) = [];

            imscheck = cell(numel(irel), 1);
            pxscheck = nan(numel(irel), 2,'single');
            for inew = 1:numel(irel)
                fprintf('Reading images from file number %d of %d. Scene %d out of %d.\n',ifile,numel(imagefiles),inew, numel(irel));
                dataim.series = irel(inew);
                sigchan  = find(strcmp(dataim.channelNames,'Cy5','Cy3','EGFP'));
                chids = [sigchan];
                imtokeep = zeros(dataim.height, dataim.width, numel(chids), 'uint16');
                for icol = 1:numel(chids)
                    imtokeep(:,:,icol) = dataim.getPlane(1, chids(icol), 1, irel(inew));
                end
                imscheck{inew} = imtokeep;
                pxscheck(inew,:) = dataim.pxSize;
            end %for each relevant scene

            for iim = 1:numel(imscheck) %resize and histogram rescale image to 8 bit then save
                currim = imscheck{iim};
                currim = single(currim);
    %             imresizefactor = unique(pxscheck(iim,:))/resizepxsize;
                saveimg = nan([ceil(size(currim,[1 2])),...
                               size(currim,3)],...
                              'single');

                for ic = 1:size(currim, 3)
                    imtosave = squeeze(currim(:, :, ic));
                    newlims  = getImageLimits(imtosave, alphas(ic));
                    imtosave = resmax*(imtosave - min(imtosave(:)))/(max(imtosave(:))-min(imtosave(:)));

                    saveimg(:,:,ic) = imtosave;
                end

                saveimg = uint8(saveimg);
                ax = subplot(1,2,1);
                imagesc(ax,squeeze(saveimg(:,:,1)))
                fsavename = sprintf('%s_series_%02d.tif', nameuse, iim);

                t = Tiff(fullfile(savefolder, fsavename),'w');
                setTag(t,'Photometric',Tiff.Photometric.MinIsBlack);
                setTag(t,'ImageLength', size(saveimg,1));
                setTag(t,'ImageWidth',  size(saveimg,2));
                setTag(t,'BitsPerSample',8);
                setTag(t,'SamplesPerPixel',1);
                setTag(t,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
                setTag(t,'ResolutionUnit',Tiff.ResolutionUnit.Centimeter);
                setTag(t,"XResolution",double(unique(pxscheck(:))*100));
                setTag(t,"YResolution",double(unique(pxscheck(:))*100));
                write(t,saveimg);
                close(t);
            end %for each relevant image rescale
            %========================================================================
            fprintf(repmat('\b', 1, numel(msg)));
            msg = sprintf('File %d/%d. Time elapsed %2.2f s...\n', ifile, numel(imagefiles),toc);
            fprintf(msg);
            %========================================================================
        end%for each file found in path
    end %if path is not empty
end %for each path

%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
%%
% for inew = 1:numel(irel)
%     dataim.series = irel(inew);
%     dataim.exportAsTIFF(finalsavepath)
% end
%%



%%
% [~, nameuse, ~] = fileparts(dpfile);
% 
% for iim = 1:numel(imscheck)
%     currim = imscheck{iim};
%     for ic = 1:size(currim, 3)
%         imtosave = single(currim(:, :, ic));
%         newlims  = getImageLimits(imtosave,0.001);
%         imtosave = 255*(imtosave - newlims(1))/range(newlims);
%         [ypix, xpix] = size(imtosave);
%         imtofit = medfilt2(imtosave, [25 25]);
%         imtofit = imtofit/max(imtofit,[],'all');
%         fitgaussrf(1:xpix, 1:ypix, double(imtofit))
% 
%         imtosave = uint8(imtosave);
%         fname = sprintf('%s_series%02d_chan%01d.png', nameuse, iim, ic);
%         imwrite(imtosave, fullfile(finalsavepath, fname))
%     end
% 
% end
% 
% 
% 
% 
% %%
% iscene = 5;
% imuv  = imscheck{iscene}(:,:,2);
% imtom = imscheck{iscene}(:,:,1);
% uvlims = getImageLimits(imuv,0.01);
% tomlim = 2*getImageLimits(imtom,0.01);
% %%
% imxlim = [2000 11500];
% imylim = [1000 7500];
% imxlimzoom = [4000 8000];
% imylimzoom = [4000 6500];
% zoomim = imtom(imylimzoom(1):imylimzoom(2), imxlimzoom(1):imxlimzoom(2));
% 
% imback    = medfilt2(zoomim, [81 81]);
% zoomimpre = single(zoomim)-single(imback);
% 
% %%
% f = figure();
% f.Units = 'centimeters';
% fw = 50; fh = 32;
% f.Position = [0 0 fw fh];
% f.MenuBar = 'none';
% f.ToolBar = 'none';
% 
% p = panel();
% p.pack('h', {0.5 0.5})
% p(1).pack('v', 2)
% p(2).pack('v', 2)
% 
% p.de.margin = 1;
% p.de.margintop = 5;
% p.fontsize = 15;
% 
% p(1,1).select();
% imagesc(imuv, uvlims);
% axis equal; xlim(imxlim); ylim(imylim);colormap(gray)
% ax = gca; ax.YDir = 'reverse'; ax.Visible = 'off';
% ax.Title.Visible = 'on';
% title('DAPI')
% p(1,2).select();
% imagesc(imtom, tomlim);colormap(gray)
% rectangle('Position',[imxlimzoom(1) imylimzoom(1) range(imxlimzoom) range(imylimzoom)],...
%     'LineWidth',2,'EdgeColor','r')
% axis equal; xlim(imxlim); ylim(imylim)
% ax = gca; ax.YDir = 'reverse'; ax.Visible = 'off';
% ax.Title.Visible = 'on';
% title('tdTomato (TRAP2 - cFOS)')
% 
% p(2,1).select();
% imagesc(zoomim,[0 quantile(zoomim,0.999,'all')]);
% axis equal; axis tight;colormap(gray)
% ax = gca; ax.YDir = 'reverse'; ax.Visible = 'off';
% ax.Title.Visible = 'on';
% title('Zoom tdTomato')
% 
% p(2,2).select();
% imagesc(zoomimpre,[0 quantile(zoomimpre,0.999,'all')]);
% axis equal; axis tight;colormap(gray)
% ax = gca; ax.YDir = 'reverse'; ax.Visible = 'off';
% ax.Title.Visible = 'on';
% title('Zoom tdTomato - local background subtraction')
% %%
% savepath = 'S:\ElboustaniLab\#SHARE\Documents\Dimos\ALiCe';
% filenamepng = sprintf('scene_%d.png',iscene);
% p.export(fullfile(savepath,filenamepng), sprintf('-w%d',fw*10),sprintf('-h%d',fh*10), '-r200');
% 
% 
% %%
% imcurr = single(imscheck{3}(:,:,1));
% imback = medfilt2(imcurr, [81 81]);
% imtop = imcurr-imback;
% 
% imuv   = single(imscheck{3}(:,:,2));
% imuv   = (imuv - min(uvlims))/range(uvlims);
% imuv(imuv>1) = 1;
% 
% imagesc(imscheck{3}(:,:,2),[0 10000]);
% colormap("bone")


% for ich = 1:2
% end
% imlims = cellfun(@(x) getImageLimits(x(:,:,1),0.01), imscheck,'UniformOutput',0);



%-------------------------------------------------------------------------
end