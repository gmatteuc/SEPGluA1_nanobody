function [outSigTiff, outKeepMaskTiff] = manual_mask(sigPath, maskPath, varargin)
% MANUAL_MASK_TIFF 
% • Draw polygons on AP (coronal) slices. 
% • AP axis (defaults to dimension == APTarget, default 900).
% • Writes NaN-masked TIFF (pages=S) and optional mask TIFF.
%
% 
% 
%     
%     

% ---------- options ----------
p = inputParser;
p.addParameter('SaveDir','',@(s)ischar(s)||isstring(s));
p.addParameter('OutName','',@(s)ischar(s)||isstring(s));
p.addParameter('Overwrite',false,@islogical);
p.addParameter('SaveMaskTiff',true,@islogical);
p.addParameter('ManualMaskDir','',@(s)ischar(s)||isstring(s));
p.addParameter('ReuseEdits',true,@islogical);
p.addParameter('StartSlice',1,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('ShowOriginalUnder',true,@islogical);
p.addParameter('Colormap','gray',@(s)ischar(s)||isstring(s));
p.addParameter('MaskThreshold',0,@(x)isnumeric(x)&&isscalar(x));
p.addParameter('MaskIndicatesKeep',true,@islogical);
p.addParameter('APFrom','auto',@(s)any(strcmpi(s,{'auto','height','pages','width'})));
p.addParameter('APTarget',900,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.parse(varargin{:});
opt = p.Results;

assert(~isempty(opt.SaveDir),'SaveDir is required.');
if ~exist(opt.SaveDir,'dir'), mkdir(opt.SaveDir); end

% ---------- signal info ----------
infoSIG = imfinfo(sigPath);
S = numel(infoSIG); H = infoSIG(1).Height; W = infoSIG(1).Width;

% ---------- load single mask or default to all-true ----------
M = read_mask_to_SHW_single(maskPath,[S H W],opt.MaskThreshold,opt.MaskIndicatesKeep);

% ---------- outputs ----------
[~,nm] = fileparts(sigPath);
baseName = char(opt.OutName); if isempty(baseName), baseName=[nm '_mask']; end
outSigTiff      = fullfile(opt.SaveDir,[baseName '_channel.tif']);
outKeepMaskTiff = fullfile(opt.SaveDir,[baseName '_binary.tif']);
if ~opt.Overwrite && exist(outSigTiff,'file') && (~opt.SaveMaskTiff || exist(outKeepMaskTiff,'file'))
    fprintf('Skipping (exists): %s\n', baseName); return;
end

% ---------- choose AP axis ----------
apFrom = choose_ap_axis(opt.APFrom, opt.APTarget, S, H, W);
switch apFrom
    case 'height', APN = H; basePlane = [S W];
    case 'pages',  APN = S; basePlane = [H W];
    case 'width',  APN = W; basePlane = [S H];
end 


polyCells = cell(APN,1);
maskDir = char(opt.ManualMaskDir);
if ~isempty(maskDir) && ~exist(maskDir,'dir'), mkdir(maskDir); end
statePath = '';
if ~isempty(maskDir)
    statePath = fullfile(maskDir,sprintf('%s_manualPolys_AP_%s.mat',baseName,apFrom));
    if opt.ReuseEdits && exist(statePath,'file')
        L = load(statePath,'polyCells');
        if isfield(L,'polyCells') && numel(L.polyCells)==APN
            polyCells = L.polyCells;
        end
    end
end

% ---------- UI ----------
cur = max(1, min(APN, round(opt.StartSlice)));
hFig = figure('Name',sprintf('Manual Masking (AP=%s)',apFrom),...
    'Color','w','Units','normalized','Position',[0.16 0.10 0.72 0.84],...
    'Visible','on','KeyPressFcn',@onKey,'WindowScrollWheelFcn',@onScroll);
ax = axes('Parent',hFig); colormap(ax,opt.Colormap); axis(ax,'image'); axis(ax,'off'); hold(ax,'on');

% Top controls
uicontrol('Style','text','String','Slice:','Units','normalized','Position',[0.02 0.94 0.05 0.04],'BackgroundColor','w');
hTxtIdx = uicontrol('Style','text','String','1/1','Units','normalized','Position',[0.07 0.94 0.10 0.04],...
    'BackgroundColor','w','HorizontalAlignment','left');
uicontrol('Style','pushbutton','String','Prev (P)','Units','normalized','Position',[0.20 0.94 0.10 0.045],...
    'Callback',@(s,e) go(-1));
uicontrol('Style','pushbutton','String','Next (N)','Units','normalized','Position',[0.31 0.94 0.10 0.045],...
    'Callback',@(s,e) go(+1));
uicontrol('Style','pushbutton','String','Add polygon (A)','Units','normalized','Position',[0.45 0.94 0.16 0.045],...
    'Callback',@(s,e) addPolygon());
uicontrol('Style','pushbutton','String','Delete selected (D)','Units','normalized','Position',[0.62 0.94 0.16 0.045],...
    'Callback',@(s,e) deleteSelected());
uicontrol('Style','pushbutton','String','Finish & Write (F)','Units','normalized','Position',[0.81 0.94 0.17 0.045],...
    'Callback',@(s,e) doneAndWrite());

% Slider
smallStep = (APN>1) * max(1/(APN-1), 1e-6);
bigStep   = (APN>1) * min(10/(APN-1), 1);
hSlider = uicontrol('Style','slider','Units','normalized','Position',[0.10 0.03 0.80 0.035],...
    'Min',1,'Max',APN,'Value',cur,'SliderStep',[smallStep bigStep],'Callback',@onSlider);
uicontrol('Style','text','String','1','Units','normalized','Position',[0.02 0.03 0.06 0.035],...
    'BackgroundColor','w','HorizontalAlignment','center');
uicontrol('Style','text','String',num2str(APN),'Units','normalized','Position',[0.91 0.03 0.06 0.035],...
    'BackgroundColor','w','HorizontalAlignment','center');

roiList = gobjects(0);
redraw(); uiwait(hFig);

    function onKey(~,evt)
        switch lower(evt.Key)
            case {'rightarrow','n'}, go(+1);
            case {'leftarrow','p'}, go(-1);
            case {'a'}, addPolygon();
            case {'d'}, deleteSelected();
            case {'f','return'}, doneAndWrite();
        end
    end
    function onScroll(~,evt)
        step = sign(evt.VerticalScrollCount);
        if step>0, go(+1); elseif step<0, go(-1); end
    end
    function onSlider(src,~)
        newVal = round(get(src,'Value'));
        if newVal ~= cur, saveCurrentPolys(); cur = clamp(newVal,1,APN); redraw();
        else, set(src,'Value',cur);
        end
    end
    function go(step)
        saveCurrentPolys();
        cur = clamp(cur+step,1,APN);
        set(hSlider,'Value',cur);
        redraw();
    end

    function redraw()
        cla(ax); axis(ax,'image','off'); hold(ax,'on');
        [img2d_base, keepBase_base] = build_base_plane(sigPath, infoSIG, M, apFrom, cur);
        if opt.ShowOriginalUnder, imagesc(ax, img2d_base);
        else, imagesc(ax, img2d_base .* single(keepBase_base)); end
        overlay(ax, ~keepBase_base, [0 1 1], 0.25); 

        roiList = gobjects(0);
        polys = polyCells{cur};
        for k = 1:numel(polys)
            P = polys{k};
            if size(P,2)==2 && size(P,1)>=3
                roi = drawpolygon(ax,'Position',P,'LineWidth',1.5,...
                    'FaceAlpha',0.25,'Color',[1 0 1],'InteractionsAllowed','all');
                roiList(end+1)=roi; 
            end
        end

        set(hTxtIdx,'String',sprintf('%d / %d',cur,APN));
        set(hSlider,'Value',cur);
        drawnow;
    end

    function addPolygon()
        roi = drawpolygon(ax,'LineWidth',1.5,'Color',[1 0 1],...
            'FaceAlpha',0.25,'InteractionsAllowed','all');
        if ~isempty(roi) && isvalid(roi)
            roiList(end+1)=roi;
        end
    end

    function deleteSelected()
        keepIdx=true(size(roiList));
        for k=1:numel(roiList)
            if isvalid(roiList(k)) && strcmp(get(roiList(k),'Selected'),'on')
                delete(roiList(k)); keepIdx(k)=false;
            end
        end
        roiList=roiList(keepIdx);
    end

    function saveCurrentPolys()
        polysOut={};
        for ii=1:numel(roiList)
            if isvalid(roiList(ii))
                P=roiList(ii).Position;
                if size(P,1)>=3, polysOut{end+1}=P; end
            end
        end
        polyCells{cur}=polysOut;
        if ~isempty(statePath), save(statePath,'polyCells','-v7.3'); end
    end

    function doneAndWrite()
        saveCurrentPolys(); try close(hFig); end
        writeOutputs(); uiresume;
    end

    function writeOutputs()
        manErase = false(S,H,W);
        switch apFrom
            case 'height'
                for h=1:H
                    BW = rasterize_polys(polyCells{h}, [S W]); % [S x W]
                    for s=1:S, manErase(s,h,:) = BW(s,:); end
                end
            case 'pages'
                for s=1:S
                    BW = rasterize_polys(polyCells{s}, [H W]); % [H x W]
                    manErase(s,:,:) = BW;
                end
            case 'width'
                for w=1:W
                    BW = rasterize_polys(polyCells{w}, [S H]); % [S x H]
                    for s=1:S, manErase(s,:,w) = BW(s,:); end
                end
        end

        % Open writers (append pages correctly without creating extra empty dirs)
        tSig=Tiff(outSigTiff,'w8'); c1=onCleanup(@()tryClose(tSig)); 
        if opt.SaveMaskTiff, tMsk=Tiff(outKeepMaskTiff,'w8'); c2=onCleanup(@()tryClose(tMsk)); else, tMsk=[]; end

        for s=1:S
            sig2d = single(imread(sigPath,s,'Info',infoSIG));
            keep  = squeeze(M(s,:,:)) & ~squeeze(manErase(s,:,:));
            sig2d(~keep) = NaN;

            % write signal slice
            if s>1, tSig.writeDirectory(); end
            setFloatTags(tSig, size(sig2d,1), size(sig2d,2));
            tSig.write(sig2d);

            % write mask slice 
            if opt.SaveMaskTiff
                if s>1, tMsk.writeDirectory(); end
                setFloatTags(tMsk, size(keep,1), size(keep,2));
                tMsk.write(single(keep));
            end
        end

        fprintf('Wrote: %s\n',outSigTiff);
        if opt.SaveMaskTiff, fprintf('Wrote: %s\n',outKeepMaskTiff); end
    end
end

% ===== helpers functions =====
function apFrom = choose_ap_axis(APFrom, APTarget, S, H, W)
if ~strcmpi(APFrom,'auto')
    apFrom = lower(APFrom); return;
end
dims = [S H W];
names = {'pages','height','width'};
idx = find(dims==APTarget);
if numel(idx)==1
    apFrom = names{idx};
else
    [~,ord] = sort(abs(dims-APTarget));   
    apFrom = names{ord(1)};               
end
end

function [img2d_base, keepBase_base] = build_base_plane(sigPath, infoSIG, M, apFrom, idx)
S = numel(infoSIG); H = infoSIG(1).Height; W = infoSIG(1).Width;
switch apFrom
    case 'height'
        h = idx; img2d_base = zeros(S,W,'single');
        for s=1:S, r=single(imread(sigPath,s,'Info',infoSIG)); img2d_base(s,:)=r(h,:); end
        keepBase_base = squeeze(M(:,h,:));
    case 'pages'
        img2d_base = single(imread(sigPath,idx,'Info',infoSIG));
        keepBase_base = squeeze(M(idx,:,:));
    case 'width'
        w = idx; img2d_base = zeros(S,H,'single');
        for s=1:S, r=single(imread(sigPath,s,'Info',infoSIG)); img2d_base(s,:)=r(:,w)'; end
        keepBase_base = squeeze(M(:,:,w));
end
end

function BW = rasterize_polys(polys, sz2d)

BW = false(sz2d);
for k = 1:numel(polys)
    P = polys{k};
    if size(P,2)==2 && size(P,1)>=3
        BW = BW | poly2mask(P(:,1), P(:,2), sz2d(1), sz2d(2));
    end
end
end

function L = read_mask_to_SHW_single(pathIn, sz, thr, maskIndicatesKeep)

S=sz(1); H=sz(2); W=sz(3);
if isempty(pathIn) || (isstring(pathIn) && strlength(pathIn)==0)
    L = true(S,H,W); return;
end
if isstring(pathIn), pathIn = char(pathIn); end
[~,~,ext]=fileparts(pathIn);
applyPol = @(X) logical(maskIndicatesKeep .* X | (~maskIndicatesKeep) .* ~X);
if any(strcmpi(ext,{'.tif','.tiff'}))
    info=imfinfo(pathIn); P=numel(info); h=info(1).Height; w=info(1).Width;
    rp=@(k) logical((imread(pathIn,k,'Info',info))>thr);
    if P==S && h==H && w==W
        L=false(S,H,W); for k=1:S, L(k,:,:) = rp(k); end; L=applyPol(L); return
    end
    if P==S && h==W && w==H
        L=false(S,H,W); for k=1:S, L(k,:,:) = rp(k)'; end
        warning('Mask %s: pages==S, H/W swapped.',pathIn); L=applyPol(L); return
    end
    if P==1 && ((h==H&&w==W)||(h==W&&w==H))
        img=rp(1); if h==W&&w==H, img=img'; end
        L=repmat(reshape(img,[1 H W]),S,1,1); L=applyPol(L); return
    end
    if P==H && h==W && w==S
        % Mask stored as [W x S] per H-slice
        L=false(S,H,W);
        for hh=1:H
            page=rp(hh); % [W x S]
            for s=1:S
                L(s,hh,:) = page(:,s)';  % ML dimension
            end
        end
        warning('Mask %s: pages==H, perpage=[W x S] → reshaped to [S H W].',pathIn);
        L=applyPol(L); return
    end
    error('Mask TIFF dims mismatch. Signal=[%d %d %d], mask pages=%d, perpage=[%d %d]',S,H,W,P,h,w);
elseif strcmpi(ext,'.mat')
    Sdata=load(pathIn); f=fieldnames(Sdata); Lraw=Sdata.(f{1});
    if ~islogical(Lraw), Lraw=Lraw>thr; end
    if isequal(size(Lraw),[S H W]), L=applyPol(Lraw); return
    end
    error('Mask MAT dims mismatch. Expected [%d %d %d].',S,H,W);
else
    error('Unsupported mask format: %s',ext);
end
end

function overlay(ax,mask,rgb,alpha)
if ~any(mask(:)), return; end
C = zeros([size(mask) 3]);
C(:,:,1)=rgb(1); C(:,:,2)=rgb(2); C(:,:,3)=rgb(3);
h = imagesc(ax,C); set(h,'AlphaData',alpha*double(mask));
end

function y = clamp(x,a,b), y = max(a, min(b, x)); end

function tryClose(tObj)
try, tObj.close(); catch, end
end

function setFloatTags(tObj, rows, cols)
tags.ImageLength      = rows;
tags.ImageWidth       = cols;
tags.Photometric      = Tiff.Photometric.MinIsBlack;
tags.Compression      = Tiff.Compression.None;
tags.SampleFormat     = Tiff.SampleFormat.IEEEFP;
tags.BitsPerSample    = 32;
tags.SamplesPerPixel  = 1;
tags.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;;
tObj.setTag(tags);
end
