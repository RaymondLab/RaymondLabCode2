frame = 434;
camera = 1;
edgeThresh = 40;
minFeatures = 0.6;
spotMaskRadius = 8;
imadj = [0.1, 0.9];
radiiPupil = [30, 80];


%% Load and process image
img = imread(fullfile(cd, ['img',num2str(camera),'.tiff']), 'Index',frame);
img = imresize(img, 2);
img = imfilter(img, fspecial('gaussian',3,.5));
img = imadjust(img, imadj, []);

%% Define estimated pupil center and roi
tempfig = figure('units','normalized', 'outerposition',[0 0 1 1]); clf
temph = tiledlayout(tempfig, 1, 1, 'TileSpacing','compact', 'Padding','compact');
tempax = nexttile;
hold(tempax, 'on');
imagesc(tempax, img);
colormap(tempax, 'gray');
axis(tempax, 'image');
set(tempax, 'ydir', 'reverse');
set(tempax, 'Visible', 'off');
center = drawcrosshair(tempax, 'LineWidth',1, 'Color','r');
roi = drawellipse(tempax, 'Color','m');
hold(tempax, 'on');

Px = center.Position(1);
Py = center.Position(2);
roimask = createMask(roi);
close(tempfig);

img(~roimask) = nan;
newimg = double(img);
imgDilated = newimg;
% Remove bright spots in image
totalMask = newimg > 204;
totalMaskDilate = imdilate(totalMask, strel('disk',spotMaskRadius));
imgDilated(totalMaskDilate) = nan;
gaussian_smooth_image = @(I,sigma) imfilter(I, fspecial('gaussian', [ceil(2.5*sigma) ceil(2.5*sigma)], sigma), 'symmetric');
imgSmoothed = gaussian_smooth_image(imgDilated, 3);

%% Detect pupil borders using starburst algorithm
[epx, epy, ~] = starburst_pupil_contour_detection(imgSmoothed, ...
                                                  Px, Py, ...
                                                  edgeThresh, ...
                                                  round(radiiPupil), ...
                                                  minFeatures);

[epxtest, epytest, ~] = starburst_pupil_contour_detection2(imgSmoothed, ...
                                                  Px, Py, ...
                                                  edgeThresh, ...
                                                  round(radiiPupil), ...
                                                  minFeatures);

[epxnew, epynew, edge_thresh_new, points.px0, points.py0] = starburstPupilContourDetection(imgSmoothed, ...
                                                                                     Px, Py, ...
                                                                                     edgeThresh, ...
                                                                                     round(radiiPupil), ...
                                                                                     minFeatures);

%% Original ellipse fitting method
[~, inliers] = fit_ellipse_ransac(epx(:), epy(:), radiiPupil + [-15, 15]);
epx2 = epx(inliers);
epy2 = epy(inliers);
ellipseResult = fit_ellipse(epx2, epy2);

pupil_x = ellipseResult.X0_in;
pupil_y = ellipseResult.Y0_in;
pupil_r1 = ellipseResult.a;
pupil_r2 = ellipseResult.b;
pupil_angle = -ellipseResult.phi;

phi = linspace(0, 2*pi, 100);
pRadius = max([pupil_r1, pupil_r2]);
pLx = pupil_r1*cos(phi)*cos(pupil_angle) - ...
        sin(pupil_angle)*pupil_r2*sin(phi) + ...
        pupil_x;
pLy = pupil_r1*cos(phi)*sin(pupil_angle) + ...
        cos(pupil_angle)*pupil_r2*sin(phi) + ...
        pupil_y;
pMx = pupil_x;
pMy = pupil_y;

%% Modified Original ellipse fitting method
[~, inlierstest] = fit_ellipse_ransac(epxtest(:), epytest(:), radiiPupil + [-15, 15]);
epx2test = epxtest(inlierstest);
epy2test = epytest(inlierstest);
ellipseResulttest = fit_ellipse(epx2test, epy2test);

pupiltest_x = ellipseResulttest.X0_in;
pupiltest_y = ellipseResulttest.Y0_in;
pupiltest_r1 = ellipseResulttest.a;
pupiltest_r2 = ellipseResulttest.b;
pupiltest_angle = -ellipseResulttest.phi;

ptestRadius = max([pupiltest_r1, pupiltest_r2]);
ptestLx = pupiltest_r1*cos(phi)*cos(pupiltest_angle) - ...
        sin(pupiltest_angle)*pupiltest_r2*sin(phi) + ...
        pupiltest_x;
ptestLy = pupiltest_r1*cos(phi)*sin(pupiltest_angle) + ...
        cos(pupiltest_angle)*pupiltest_r2*sin(phi) + ...
        pupiltest_y;
ptestMx = pupiltest_x;
ptestMy = pupiltest_y;

% %% Hannah's updated ellipse fitting method
% p_robustEllipse = [];
% p_robustEllipse.min_dist = 2; % distance threshold for segmentation (t)
% p_robustEllipse.min_num = 2;  % min number of points per cluster
% p_robustEllipse.D = 10;       % expected number of sets
% p_robustEllipse.S = 3;        % number of additional sets to keep in each iteration
% p_robustEllipse.S_max = 10;   % max number of subsets to keep in each iteration
% p_robustEllipse.sigma = 1.15; % error threshold; searching process converges if excluding any subset does not reduce the energy up to certain rate σ
% p_robustEllipse.eta = 2.5;      % default 5. error threshold; in case there are much more outliers than average subset size so that removing any subset does not reduce the energy
% p_robustEllipse.inclusion_penalty = 6;
% p_robustEllipse.size_penalty = 10;
% p_robustEllipse.size_ellipse = mean(radiiPupil); % Guess to enforce correct pupil size (pixels)
% [points.px1, points.py1, ellipseResult2] = robustEllipse(epxnew, epynew, p_robustEllipse);
% 
% pupil2_x = ellipseResult2.x0;
% pupil2_y = ellipseResult2.y0;
% pupil2_r1 = ellipseResult2.a;
% pupil2_r2 = ellipseResult2.b;
% pupil2_angle = -ellipseResult2.phi;
% 
% p2Radius = max([pupil2_r1, pupil2_r2]);
% p2Lx = pupil2_r1*cos(phi)*cos(pupil2_angle) - ...
%         sin(pupil2_angle)*pupil2_r2*sin(phi) + ...
%         pupil2_x;
% p2Ly = pupil2_r1*cos(phi)*sin(pupil2_angle) + ...
%         cos(pupil2_angle)*pupil2_r2*sin(phi) + ...
%         pupil2_y;
% p2Mx = pupil2_x;
% p2My = pupil2_y;

% %% Do better fit of resulting points
% ellipseResult3 = fitEllipse(points.px1, points.py1);
% 
% pupil3_x = ellipseResult3.x0;
% pupil3_y = ellipseResult3.y0;
% pupil3_r1 = ellipseResult3.a;
% pupil3_r2 = ellipseResult3.b;
% pupil3_angle = -ellipseResult3.phi;
% 
% p3Radius = max([pupil3_r1, pupil3_r2]);
% p3Lx = pupil3_r1*cos(phi)*cos(pupil3_angle) - ...
%         sin(pupil3_angle)*pupil3_r2*sin(phi) + ...
%         pupil3_x;
% p3Ly = pupil3_r1*cos(phi)*sin(pupil3_angle) + ...
%         cos(pupil3_angle)*pupil3_r2*sin(phi) + ...
%         pupil3_y;
% p3Mx = pupil3_x;
% p3My = pupil3_y;
% 
% 
% %% Original ellipse fitting method with new starburst function
% [~, inliersnew] = fit_ellipse_ransac(epxnew(:), epynew(:), radiiPupil + [-15, 15]);
% epx2new = epxnew(inliersnew);
% epy2new = epynew(inliersnew);
% ellipseResultnew = fit_ellipse(epx2new, epy2new);
% 
% pupilnew_x = ellipseResultnew.X0_in;
% pupilnew_y = ellipseResultnew.Y0_in;
% pupilnew_r1 = ellipseResultnew.a;
% pupilnew_r2 = ellipseResultnew.b;
% pupilnew_angle = -ellipseResultnew.phi;
% 
% pnewRadius = max([pupilnew_r1, pupilnew_r2]);
% pnewLx = pupilnew_r1*cos(phi)*cos(pupilnew_angle) - ...
%         sin(pupilnew_angle)*pupilnew_r2*sin(phi) + ...
%         pupilnew_x;
% pnewLy = pupilnew_r1*cos(phi)*sin(pupilnew_angle) + ...
%         cos(pupilnew_angle)*pupilnew_r2*sin(phi) + ...
%         pupilnew_y;
% pnewMx = pupilnew_x;
% pnewMy = pupilnew_y;


%% Plot
fig = figure('units','normalized', 'outerposition',[0 0 1 1]); clf
h = tiledlayout(fig, 1, 2, 'TileSpacing','compact', 'Padding','compact');

ax1 = nexttile;
imagesc(ax1, imgSmoothed);
hold(ax1, 'on');
line(ax1, pLx, pLy, 'Color','r', 'LineWidth',2);
plot(ax1, pMx, pMy, '+r', 'LineWidth',2, 'MarkerSize',10);
plot(ax1, epx, epy, '.y');
plot(ax1, epx2, epy2, '.c');
colormap(ax1, 'gray');
axis(ax1, 'image');
set(ax1, 'ydir','reverse');
set(ax1, 'Visible','off');
title(ax1, "Original Method: OLD Starburst + fit_ellipse_ransac + fit_ellipse", 'Interpreter','none');
ax1.Title.Visible = 'on';
hold(ax1, 'off');

ax2 = nexttile;
imagesc(ax2, imgSmoothed);
hold(ax2, 'on');
line(ax2, ptestLx, ptestLy, 'Color','r', 'LineWidth',2);
plot(ax2, ptestMx, ptestMy, '+r', 'LineWidth',2, 'MarkerSize',10);
plot(ax2, epxtest, epytest, '.y');
plot(ax2, epx2test, epy2test, '.c');
colormap(ax2, 'gray');
axis(ax2, 'image');
set(ax2, 'ydir','reverse');
set(ax2, 'Visible','off');
title(ax2, 'Test Method 1: Modified Old Starburst + fit_ellipse_ransac + fit_ellipse', 'Interpreter','none');
ax2.Title.Visible = 'on';
hold(ax2, 'off');

% ax3 = nexttile;
% imagesc(ax3, imgSmoothed);
% hold(ax3, 'on');
% line(ax3, p3Lx, p3Ly, 'Color','r', 'LineWidth',2);
% plot(ax3, p3Mx, p3My, '+r', 'LineWidth',2, 'MarkerSize',10);
% plot(ax3, epxnew, epynew, '.y');
% plot(ax3, points.px1, points.py1, '.c');
% colormap(ax3, 'gray');
% axis(ax3, 'image');
% set(ax3, 'ydir','reverse');
% set(ax3, 'Visible','off');
% title(ax3, 'New Method 2: NEW Starburst + robustEllipse + fitEllipse', 'Interpreter','none');
% ax3.Title.Visible = 'on';
% hold(ax3, 'off');
% 
% ax4 = nexttile;
% imagesc(ax4, imgSmoothed);
% hold(ax4, 'on');
% line(ax4, pnewLx, pnewLy, 'Color','r', 'LineWidth',2);
% plot(ax4, pnewMx, pnewMy, '+r', 'LineWidth',2, 'MarkerSize',10);
% plot(ax4, epxnew, epynew, '.y');
% plot(ax4, epx2new, epy2new, '.c');
% colormap(ax4, 'gray');
% axis(ax4, 'image');
% set(ax4, 'ydir','reverse');
% set(ax4, 'Visible','off');
% title(ax4, "New + Old Method: NEW Starburst + fit_ellipse_ransac + fit_ellipse", 'Interpreter','none');
% ax4.Title.Visible = 'on';
% hold(ax4, 'off');

linkaxes([ax1, ax2], 'xy'); %, ax3, ax4