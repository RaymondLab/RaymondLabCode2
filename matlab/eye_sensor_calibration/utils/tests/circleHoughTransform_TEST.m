% Camera and frame number of image to load
cam = 1;
frame = 1;

% imageAdjust low/high parameters
lowIn = 0.1; %0.1;
highIn = 0.85; %0.85;

% Colorbar limit
climit = 255;


%% Script starts here

% Load sample image
img = imread(fullfile(cd,['img',num2str(cam),'.tiff']),'Index',frame);

% Upsample image
img = imresize(img, 2);

% Contrast adjusted and normalized cubic transformed images
imgAdjust = imadjust(img,[lowIn,highIn],[]);
imgCRs = double(imgAdjust);

[A, circen, crr] = CircularHough_Grd_TEST(imgCRs, [6, 13], 10, 8, 1);

%% Plot figure
fig = figure('units','normalized', 'outerposition',[0.05 0.2 0.9 0.6]);

ax(1) = subplot(1, 2, 1, 'Parent',fig);
ax1 = imagesc(ax(1), imgCRs);
colormap(ax(1), 'gray');
set(ax(1), 'ydir', 'reverse');
set(ax(1), 'Visible', 'off');

ax(2) = subplot(1, 2, 2, 'Parent',fig);
ax2 = imagesc(ax(2), A);
hold(ax(2), 'on');
for i = 1:length(crr)
    plot(ax(2), circen(i,1), circen(i,2), 'rx', 'LineWidth',3, 'MarkerSize',12);
    plotcircle(ax(2), circen(i,1), circen(i,2), crr(i), 'r');
end
colormap(ax(2), 'parula');
colorbar(ax(2));
set(ax(2), 'ydir', 'reverse');
set(ax(2), 'Visible', 'off');

linkaxes([ax(1), ax(2)], 'xy');



function h = plotcircle(target, x, y, r, color)
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(target, xunit, yunit, color, 'LineWidth',2);
end

function dispMinMaxText(name, image)
    minmaxtext = sprintf([name,' minmax = (',num2str(min(image(:))),', ', num2str(max(image(:))),')']);
    disp(minmaxtext);
end