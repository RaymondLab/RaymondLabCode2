function saveGroundTruthFrame(vars, trackParams, frame)
%SAVEGROUNDTRUTHFRAME

[~,filename] = fileparts(vars.currentFolderPath);

% Get resolution of smallest available display
monitors = get(groot,'MonitorPositions');
[~, minIdx] = min(prod(monitors(:,3:end), 2));
maxwh = monitors(minIdx, 3:4);
figwh = round(maxwh.*0.7);
figxy = round((maxwh-figwh).*0.5);
uixywh = [figxy, figwh];
clear monitors minIdx figwh figxy

% Initialize main window
fig = uifigure('Toolbar','none', ... 
               'MenuBar','none', ...
               'NumberTitle','off', ...
               'Name','Save Ground Truth Frame', ...
               'Position',uixywh, ...
               'WindowState','maximized');

% Define main window grid
g0 = uigridlayout(fig, ...
    'RowHeight',{100, '1x', 90}, ...
    'ColumnWidth',{'1x'});

% Set up panel for shape tools
pbuttons1 = uipanel(g0);
gpbuttons1 = uigridlayout(pbuttons1, ...
    'RowHeight',{'1x'}, ...
    'ColumnWidth',{200, 200, 200, '1x'}, ...
    'ColumnSpacing',20, ...
    'Padding',[20 20 20 20]);
bpupil = uibutton(gpbuttons1, ...
    'state', ...
    'Text','PUPIL');
bcr1 = uibutton(gpbuttons1, ...
    'state', ...
    'Text','CR 1');
bcr2 = uibutton(gpbuttons1, ...
    'state', ...
    'Text','CR 2');

% Set up panel for displayed frames
pplots = uipanel(g0);
gpplots = uigridlayout(pplots, ...
    'RowHeight',{'1x'}, ...
    'ColumnWidth',{'1x', '1x', '1x'}, ...
    'ColumnSpacing',20, ...
    'Padding',[20 20 20 20]);
ax1 = uiaxes(gpplots);
ax2 = uiaxes(gpplots);
ax3 = uiaxes(gpplots);

% Set up panel for save and cancel buttons
pbuttons2 = uipanel(g0);
gpbuttons2 = uigridlayout(pbuttons2, ...
    'RowHeight',{'1x'}, ...
    'ColumnWidth',{200, 200, '1x'}, ...
    'ColumnSpacing',20, ...
    'Padding',[20 20 20 20]);
bsave = uibutton(gpbuttons2, ...
    'Text','Save');
bcancel = uibutton(gpbuttons2, ...
    'Text','Cancel');

% Preallocate shape data struct
shapes = cell(1, 3);
imgdata = cell(1, 3);
for ii = 1:3
    shapes{ii} = struct;
    shapes{ii}.pupilEllipse = [];
    shapes{ii}.cr1Circle = [];
    shapes{ii}.cr2Circle = [];

    imgdata{ii} = struct;
    imgdata{ii}.pupil_x = nan;
    imgdata{ii}.pupil_y = nan;
    imgdata{ii}.pupil_r1 = nan;
    imgdata{ii}.pupil_r2 = nan;
    imgdata{ii}.pupil_angle = nan;

    imgdata{ii}.cr1_x = nan;
    imgdata{ii}.cr1_y = nan;
    imgdata{ii}.cr1_r = nan;

    imgdata{ii}.cr2_x = nan;
    imgdata{ii}.cr2_y = nan;
    imgdata{ii}.cr2_r = nan;
end

cam = trackParams.cam;
minFrames = 1;
maxFrames = round(trackParams.nImages - 1);

img = imread(fullfile(cd, ['img',num2str(cam),'.tiff']), 'Index',frame);
img = imresize(img, 2);
imgs = zeros([size(img) 3], 'uint8');
imgs(:,:,2) = img;
cxy = round(size(img).*0.5);

% Set up left frame image
if (frame-1) >= minFrames
    img = imread(fullfile(cd, ['img',num2str(cam),'.tiff']), 'Index',frame-1);
    img = imresize(img, 2);
    imgs(:,:,1) = img;
    imagesc(ax1, imgs(:,:,1));
end
title(ax1, sprintf(['Frame ',num2str(frame-1)]), 'FontSize',18);
colormap(ax1, 'gray');
axis(ax1, 'image');
set(ax1, 'ydir', 'reverse');
set(ax1, 'Visible', 'off');
set(get(ax1,'Title'), 'Visible','on');

% Set up middle (main) frame image
imagesc(ax2, imgs(:,:,2));
title(ax2, sprintf(['Frame ',num2str(frame)]), 'FontSize',18);
colormap(ax2, 'gray');
axis(ax2, 'image');
set(ax2, 'ydir', 'reverse');
set(ax2, 'Visible', 'off');
set(get(ax2,'Title'), 'Visible','on');

% Set up left frame image
if (frame-1) <= maxFrames
    img = imread(fullfile(cd, ['img',num2str(cam),'.tiff']), 'Index',frame+1);
    img = imresize(img, 2);
    imgs(:,:,3) = img;
    imagesc(ax3, imgs(:,:,3));
end
title(ax3, sprintf(['Frame ',num2str(frame+1)]), 'FontSize',18);
colormap(ax3, 'gray');
axis(ax3, 'image');
set(ax3, 'ydir', 'reverse');
set(ax3, 'Visible', 'off');
set(get(ax3,'Title'), 'Visible','on');

% Set interactivity functions for buttons
bpupil.ValueChangedFcn  = @(src,evt) bpupilPushed(evt, cxy, ax1, ax2, ax3);
bcr1.ValueChangedFcn    = @(src,evt) bcr1Pushed(evt, cxy, ax1, ax2, ax3);
bcr2.ValueChangedFcn    = @(src,evt) bcr2Pushed(evt, cxy, ax1, ax2, ax3);
bsave.ButtonPushedFcn   = @(src,evt) bsavePushed(fig, filename, frame, imgs);
bcancel.ButtonPushedFcn = @(src,evt) bcancelPushed(fig);

% Wait until figure is closed to exit function
uiwait(fig);

    function bpupilPushed(evt, cxy, ax1, ax2, ax3)
    if evt.Value == 1 && isempty(shapes{1}.pupilEllipse)
        shapes{1}.pupilEllipse = drawellipse(ax1, 'Center',cxy, 'SemiAxes',[75,50], 'Color','r', 'LineWidth',1);
        shapes{2}.pupilEllipse = drawellipse(ax2, 'Center',cxy, 'SemiAxes',[75,50], 'Color','r', 'LineWidth',1);
        shapes{3}.pupilEllipse = drawellipse(ax3, 'Center',cxy, 'SemiAxes',[75,50], 'Color','r', 'LineWidth',1);
    elseif evt.Value == 1 
        for nn = 1:3
            shapes{nn}.pupilEllipse.Visible = 'on';
        end
    else
        for nn = 1:3
            shapes{nn}.pupilEllipse.Visible = 'off';
        end
    end
    end

    function bcr1Pushed(evt, cxy, ax1, ax2, ax3)
    if evt.Value == 1 && isempty(shapes{1}.cr1Circle)
        shapes{1}.cr1Circle = drawcircle(ax1, 'Center',cxy, 'Radius',50, 'Color','b');
        shapes{2}.cr1Circle = drawcircle(ax2, 'Center',cxy, 'Radius',50, 'Color','b');
        shapes{3}.cr1Circle = drawcircle(ax3, 'Center',cxy, 'Radius',50, 'Color','b');
    elseif evt.Value == 1
        for nn = 1:3
            shapes{nn}.cr1Circle.Visible = 'on';
        end
    else
        for nn = 1:3
            shapes{nn}.cr1Circle.Visible = 'off';
        end
    end
    end

    function bcr2Pushed(evt, cxy, ax1, ax2, ax3)
    if evt.Value == 1 && isempty(shapes{1}.cr2Circle)
        shapes{1}.cr2Circle = drawcircle(ax1, 'Center',cxy, 'Radius',50, 'Color','c');
        shapes{2}.cr2Circle = drawcircle(ax2, 'Center',cxy, 'Radius',50, 'Color','c');
        shapes{3}.cr2Circle = drawcircle(ax3, 'Center',cxy, 'Radius',50, 'Color','c');
    elseif evt.Value == 1
        for nn = 1:3
            shapes{nn}.cr2Circle.Visible = 'on';
        end
    else
        for nn = 1:3
            shapes{nn}.cr2Circle.Visible = 'off';
        end
    end
    end

    function bsavePushed(fig, filename, frame, imgs)
    nopupil = any(cellfun(@(x) isempty(x.pupilEllipse), shapes));
    nocr1 = any(cellfun(@(x) isempty(x.cr1Circle), shapes));
    nocr2 = any(cellfun(@(x) isempty(x.cr2Circle), shapes));
    if any([nopupil, nocr1, nocr2])
        return
    end
    for nn = 1:3
        imgdata{nn}.pupil_x = shapes{nn}.pupilEllipse.Center(1);
        imgdata{nn}.pupil_y = shapes{nn}.pupilEllipse.Center(2);
        imgdata{nn}.pupil_r1 = shapes{nn}.pupilEllipse.SemiAxes(2);
        imgdata{nn}.pupil_r2 = shapes{nn}.pupilEllipse.SemiAxes(1);
        imgdata{nn}.pupil_angle = shapes{nn}.pupilEllipse.RotationAngle;

        imgdata{nn}.cr1_x = shapes{nn}.cr1Circle.Center(1);
        imgdata{nn}.cr1_y = shapes{nn}.cr1Circle.Center(2);
        imgdata{nn}.cr1_r = shapes{nn}.cr1Circle.Radius;
    
        imgdata{nn}.cr2_x = shapes{nn}.cr2Circle.Center(1);
        imgdata{nn}.cr2_y = shapes{nn}.cr2Circle.Center(2);
        imgdata{nn}.cr2_r = shapes{nn}.cr2Circle.Radius;
    end
    newfilename = sprintf('%s_%s.mat', filename, num2str(frame));
    uisave({'imgs','imgdata'}, newfilename);
    close(fig);
    end

    function bcancelPushed(fig)
    close(fig);
    end

end