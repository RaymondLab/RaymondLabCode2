%{
Relevant files in Hannah's script folder:

processEyetrack.m : Has initial parameters
trackEye.m        : Has call to next script
detectPupilCR.m   : Similar to our version of the script
%}

clear;

%% Parameters
plot_1 = true;
plot_2 = true;

% Raymond Lab parameters
frame = 1030;
camera = 1;
edgeThresh = 40;
minFeatures = 0.6;
spotMaskRadius = 8;
imadj = [0.15, 0.75];
radiiPupil = [30, 60];
alphaFRST = 2;  % Default: 0.5;
edgePadPct = 0.1;

% Hannah's parameters
hp.debug_on        = 0;
hp.smoothSigma     = 3;      % (pixels) smooth image (Default: 3)
hp.edge_thresh0    = 128;    % Initial guess, increase if pupil fit is too small, decrease to run faster (Default: 128)
hp.pupil_intensity = 20;     % Adjust based on images if needed % 30 % 100/130 (Default: 20)
hp.pupil_downsamp  = 4;      % (pixels) factor to downsample for initial pupil search (Default: 4)
hp.radiiPupil      = 55;     % (pixels) Controls approximate radius of radial symmetry filter (Default: 32)
hp.pupil_alpha     = 2;      % Sensitivity to radial symmetry. 1 - slack, 3 - strict, ... and you can go higher (Default: 2)
hp.CR_box_width    = 100;    % (pixels) Crop image for faster CR search (Default: 100)
hp.radiiCR         = [2 5];  % (pixels) [4 7] (Default: [2 5])
hp.CRthresh        = 12;     % Threshold for subtracting background (Default: 12)
hp.CRfilter        = 6;      % (default 3, minimum 3). Set larger to detect imperfect circles - but sometimes too large throws an error esp with small radiiCR
hp.minfeatures     = 0.6;     % (fraction, [0 1]) min features required for starburst pupil detection (Default: 0.9)
hp.ellipse_residual_thresh = 1.3;  % Max mean residual to accept for final ellipse fit (Default: 1.3)
hp.nCRs = 2;

p_robustEllipse         = [];
p_robustEllipse.min_num = 2;     % min number of points per cluster
p_robustEllipse.sigma   = 1.15;  % error threshold; searching process converges if excluding any subset does not reduce the energy up to certain rate σ
p_robustEllipse.eta     = 2.5;    % default 5. error threshold; in case there are much more outliers than average subset size so that removing any subset does not reduce the energy
p_robustEllipse.inclusion_penalty = 6;
p_robustEllipse.size_penalty = 10;
p_robustEllipse.size_ellipse = mean(hp.radiiPupil); % Guess to enforce correct pupil size (pixels)



%% Load and process image
imgUINT8 = imread(fullfile(cd, ['img',num2str(camera),'.tiff']), 'Index',frame);
imgUINT8 = imresize(imgUINT8, 2);
imgSINGLE = single(imgUINT8);



%% Our pipeline for Pupil Centers: 
imgRL = imfilter(imgUINT8, fspecial('gaussian',3,.5));
imgRL = imadjust(imgRL, imadj, []);

% Compute radial symmetry transform
imgRL_pupilRadial = Radial_Sym_Transform(imgRL, radiiPupil, alphaFRST);

% Crop away edges (which tend to be dark due to edge artifacts)
imgPad = round(size(imgRL_pupilRadial).*edgePadPct);
roiMask = true(size(imgRL));
roiMask(1:imgPad(1),:) = false;
roiMask(end-imgPad(1):end,:) = false;
roiMask(:,1:imgPad(2)) = false;
roiMask(:,end-imgPad(2):end) = false;
imgRL_pupilRadial(~roiMask) = nan;

% Find the minimum for pupil centers
[pupilY(1), pupilX(1)] = find(min(imgRL_pupilRadial(:))==imgRL_pupilRadial);



%% Hannah's pipeline for Pupil Centers: 
imgHP = imgaussfilt(imgSINGLE, hp.smoothSigma);
imgHP_cr = imgHP - hp.pupil_intensity;
imgHP_cr = imgHP_cr/152*255;  % Temporary hard coding
imgHP_pupil                  = imgHP_cr*3;  % For pupil detection, brighten image further so iris ~255*3/4 % *3
imgHP_pupil(imgHP_pupil>255) = 255;

hp.border_remove = 0.05;
hp.pupil_start = [NaN NaN];
imgHP_pupilCrop = imgHP_pupil(1:hp.pupil_downsamp:end, 1:hp.pupil_downsamp:end);
imgHP_pupilRadial = radialSymTransform(imgHP_pupilCrop, round(hp.radiiPupil/hp.pupil_downsamp * [.9 1.1]), hp.pupil_alpha);  
% alpha - radial strictness parameter.
%       1 - slack, accepts features with bilateral symmetry.
%       2 - a reasonable compromise.
%       3 - strict, only accepts radial symmetry.

imgHP_pupilRadial = removeBorder(imgHP_pupilRadial, hp.border_remove);
[minHP_radialPupil, indHP] = min(imgHP_pupilRadial(:));
[pupilY(2), pupilX(2)] = ind2sub(size(imgHP_pupilCrop), indHP);
pupilY(2) = pupilY(2)*hp.pupil_downsamp;
pupilX(2) = pupilX(2)*hp.pupil_downsamp;



%% Plot Pupil Center Results
if plot_1
    fig1 = figure('units','normalized', 'outerposition',[0 0 1 1]); clf
    h1 = tiledlayout(fig1, 2, 3, 'TileSpacing','tight', 'Padding','compact');
    
    plotImgs = {imgUINT8, imgRL, imgRL_pupilRadial, ...
                imgHP, imgHP_pupil, imgHP_pupilRadial};
    colLabels = ["General preprocessing", "Pupil Processing", "Radial Transformed"];
    rowIdx = 0;
    for i = 1:6
        axs1(i) = axes(h1, 'XTick',[], 'YTick',[]);
        axs1(i).Layout.Tile = i;
        hold(axs1(i), 'on');
        imagesc(axs1(i), plotImgs{i});
        if mod(i, 3) == 0
            colormap(axs1(i), 'parula');
        else
            colormap(axs1(i), 'gray');
        end
        axis(axs1(i), 'image');
        set(axs1(i), 'ydir', 'reverse');
        if mod(i, 3) == 1
            rowIdx = rowIdx + 1;
            plot(axs1(i), pupilX(rowIdx), pupilY(rowIdx), '+r', 'MarkerSize',24, 'LineWidth',3);
        end
        if rowIdx == 1
            title(axs1(i), colLabels(i), 'FontSize',14);
        end
        hold(axs1(i), 'off');
    end
    
    ylabel(axs1(1), 'Raymond Lab Version', 'FontSize',14);
    ylabel(axs1(4), 'Hannah Payne Version', 'FontSize',14);
    title(h1, 'Comparison of Pupil Center Detection', 'FontSize',18);
end



%% Our pipeline for Pupil Ellipse Fitting: 
imgRL_pupil = imgRL;
imgRL_pupil(~roiMask) = nan;
imgRL_double = double(imgRL_pupil);
imgRL_dilated = imgRL_double;

% Remove bright spots in image
totalMask = imgRL_double > 204;
totalMaskDilate = imdilate(totalMask, strel('disk',spotMaskRadius));
imgRL_dilated(totalMaskDilate) = nan;
gaussian_smooth_image = @(I,sigma) imfilter(I, fspecial('gaussian', [ceil(2.5*sigma) ceil(2.5*sigma)], sigma), 'symmetric');
imgRL_smoothed = gaussian_smooth_image(imgRL_dilated, 3);

imgRadialDiff = cell(1, 2);
imgRadialDiff{1} = computeRadialDifferences(imgRL_smoothed, pupilX(1), pupilY(1));

% Detect pupil borders using starburst algorithm
epx = cell(1, 2);
epy = cell(1, 2);
[epx{1}, epy{1}, ~] = starburst_pupil_contour_detection(imgRL_smoothed, ...
                                                        pupilX(1), pupilY(1), ...
                                                        edgeThresh, ...
                                                        round(radiiPupil), ...
                                                        minFeatures);
% Ellipse fitting method
epx2 = cell(1, 2);
epy2 = cell(1, 2);
inliers = cell(1, 2);
ellipseResult = cell(1, 2);
pupilData = cell(1, 2);

[~, inliers{1}] = fit_ellipse_ransac(epx{1}(:), epy{1}(:), radiiPupil + [-15, 15]);
epx2{1} = epx{1}(inliers{1});
epy2{1} = epy{1}(inliers{1});
ellipseResult{1} = fit_ellipse(epx2{1}, epy2{1});

pupilData{1} = getPupilData(ellipseResult{1}.X0_in, ...
                            ellipseResult{1}.Y0_in, ...
                            ellipseResult{1}.a, ...
                            ellipseResult{1}.b, ...
                            -ellipseResult{1}.phi);



%% Hannah's pipeline for Pupil Ellipse Fitting: 

% Find corneal reflections
r = hp.CR_box_width;
start_x = max(1,pupilX(2)-r/2);
stop_x = min(size(imgHP_cr, 2), pupilX(2) + r/2);
start_y = max(1,pupilY(2)-r/2); % rows = y
stop_y = min(size(imgHP_cr, 1), pupilY(2) + r/2);
imgHP_crCrop  = imgHP_cr(start_y:stop_y, start_x:stop_x); % cols = x
% Slower but better:
[~, crxy0, crr0] = CircularHoughGrd(imgHP_crCrop, hp.radiiCR, hp.CRthresh, hp.CRfilter, 1);
crxy0 = bsxfun(@plus, crxy0, [start_x start_y]-1);

% Remove CR glints
imgHP_mask = false(size(imgHP_pupil));
imgHP_mask(round(crxy0(:,2)), round(crxy0(:,1))) = true;
imgHP_mask = imdilate(imgHP_mask,strel('disk', round(mean(hp.radiiCR)*1.3))); % Expand the mask
imgHP_pupilNoCR = imgHP_pupil;
imgHP_pupilNoCR(imgHP_mask) = NaN;
imgHP_pupil = imgHP_pupilNoCR;

imgRadialDiff{2} = computeRadialDifferences(imgHP_pupil, pupilX(2), pupilY(2));

% Find pupil
mean_residual = NaN;
[epx{2}, epy{2}, edge_thresh_new, points.px0, points.py0] = starburstPupilContourDetection(imgHP_pupil, ...
    pupilX(2), pupilY(2), hp.edge_thresh0, round(hp.radiiPupil * [.8 1.2]), hp.minfeatures);

if length(points.px0) > 20  % Brian: changed points.px0 to epx{2}
    % [points.px1, points.py1, ellipse_result1] = robustEllipse(points.px0, points.py0, p_robustEllipse, hp.debug_on);
    % 
    % % Do better fit of resulting points
    % ellipseResult{2} = fitEllipse(points.px1, points.py1);
    % 
    % pupilData{2} = getPupilData(ellipseResult{2}.x0, ...
    %                             ellipseResult{2}.y0, ...
    %                             ellipseResult{2}.a, ...
    %                             ellipseResult{2}.b, ...
    %                             -ellipseResult{2}.phi);
    % epx2{2} = points.px1';
    % epy2{2} = points.py1';

    % [~, inliers{2}] = fit_ellipse_ransac(epx{2}(:), epy{2}(:), radiiPupil + [-15, 15]);
    [~, inliers{2}] = fit_ellipse_ransac(points.px0(:), points.py0(:), radiiPupil + [-15, 15]);
    epx2{2} = epx{2}(inliers{2});
    epy2{2} = epy{2}(inliers{2});
    ellipseResult{2} = fit_ellipse(epx2{1}, epy2{1});
    pupilData{2} = getPupilData(ellipseResult{2}.X0_in, ...
                                ellipseResult{2}.Y0_in, ...
                                ellipseResult{2}.a, ...
                                ellipseResult{2}.b, ...
                                -ellipseResult{2}.phi);

    % % Find the ellipse residuals
    % a = linspace(0,2*pi,200);
    % ellipse_x = ellipseResult{2}.a*cos(a)*cos( ellipseResult{2}.phi) - sin( ellipseResult{2}.phi)*ellipseResult{2}.b*sin(a) + ellipseResult{2}.x0;
    % ellipse_y = ellipseResult{2}.a*cos(a)*sin( ellipseResult{2}.phi) + cos( ellipseResult{2}.phi)*ellipseResult{2}.b*sin(a) + ellipseResult{2}.y0;
    % residuals = NaN(length(points.px1),1);
    % for ii = 1:length(points.px1)
    %     dist_sq = ((points.px1(ii) - ellipse_x(:)).^2 + (points.py1(ii) - ellipse_y(:)).^2);
    %     residuals(ii) = sqrt(min(dist_sq));
    % end
    % mean_residual = mean(residuals);
    % if mean_residual < hp.ellipse_residual_thresh
    %     pupilHP = [ellipseResult{2}.x0 ellipseResult{2}.y0 ellipseResult{2}.a ellipseResult{2}.b ellipseResult{2}.phi];
    % end
    % 
    % % Calculate the fraction of original point retained in the final fit
    % points_fraction = length(points.px0)/length(points.px1);
end

% % Check for CR duplicates
% i = 1;
% crxy = [];
% crr = [];
% crr0_init = crr0;
% crxy0_init = crxy0;
% while ~isempty(crr0_init)%i<length(crr)
%     distThresh = 10; % (px) minimum separation between CRs
%     duplicates = sqrt(sum((repmat(crxy0_init(i,:),size(crxy0_init,1),1) - crxy0_init).^2,2)) < distThresh; % 30
%     crxy(end+1,:) = mean(crxy0_init(duplicates,:),1);
%     crr(end+1,:) = mean(crr0_init(duplicates,:),1);
%     crr0_init(duplicates) = [];
%     crxy0_init(duplicates,:) = [];
% end
% 
% % If more than N CRs detected, sort by distance from pupil
% if size(crxy,1) > hp.nCRs
%     dist = sqrt(sum((bsxfun(@minus, crxy, pupilHP(1:2))).^2, 2));
%     [~, inds] = sort(dist);
%     crxy = crxy(inds(1:hp.nCRs),:);
%     crr = crr(inds(1:hp.nCRs),:);
% end
% 
% % Sort remaining CRs by horizontal position
% [crx, ind] = sort(crxy(:,1),1,'ascend');
% cry = crxy(ind,2);
% crr = crr(ind);
% 
% % Fill in other CRs with NaNs if not present
% CRx(1:length(crx)) = crx;
% CRy(1:length(cry)) = cry;
% CRr(1:length(crr)) = crr;
% CRout = [CRx; CRy; CRr];



%% Plot Pupil Detection Results
if plot_2
    fig2 = figure('units','normalized', 'outerposition',[0 0 1 1]); clf
    h2 = tiledlayout(fig2, 2, 2, 'TileSpacing','compact', 'Padding','compact');
    
    plotImgs = {imgRL_double, imgRadialDiff{1}, ...
                imgHP, imgRadialDiff{2}};
    colLabels = ["Initial image", "Detected pupil"];
    nColsPerRow = length(colLabels);
    rowIdx = 0;
    for i = 1:round(2*nColsPerRow)
        axs2(i) = axes(h2, 'XTick',[], 'YTick',[]);
        axs2(i).Layout.Tile = i;
        hold(axs2(i), 'on');
        imagesc(axs2(i), plotImgs{i});
        if mod(i, nColsPerRow) == 0
            colormap(axs2(i), 'parula');
        else
            colormap(axs2(i), 'gray');
        end
        axis(axs2(i), 'image');
        set(axs2(i), 'ydir', 'reverse');
        if mod(i, nColsPerRow) == 1
            rowIdx = rowIdx + 1;
            line(axs2(i), pupilData{rowIdx}.Lx, pupilData{rowIdx}.Ly, 'Color','r', 'LineWidth',2);
            plot(axs2(i), pupilData{rowIdx}.x, pupilData{rowIdx}.y, '+r', 'LineWidth',2, 'MarkerSize',14);
            plot(axs2(i), epx{rowIdx}, epy{rowIdx}, '.y');
            plot(axs2(i), epx2{rowIdx}, epy2{rowIdx}, '.c', 'MarkerSize',10);
        end
        if rowIdx == 1
            title(axs2(i), colLabels(i), 'FontSize',14);
        end
        hold(axs2(i), 'off');
    end
    
    ylabel(axs2(1), 'Raymond Lab Version', 'FontSize',14);
    ylabel(axs2(3), 'Hannah Payne Version', 'FontSize',14);
    title(h2, 'Comparison of Pupil Edge Detection', 'FontSize',18);
end
    
function pupil = getPupilData(x, y, r1, r2, angle)
    pupil.x = x;
    pupil.y = y;
    pupil.r1 = r1;
    pupil.r2 = r2;
    pupil.angle = angle;
    
    phi = linspace(0, 2*pi, 100);
    pupil.radius = max([r1, r2]);
    pupil.Lx = r1*cos(phi)*cos(angle) - sin(angle)*r2*sin(phi) + x;
    pupil.Ly = r1*cos(phi)*sin(angle) + cos(angle)*r2*sin(phi) + y;
end



% [epxtest, epytest, ~] = starburst_pupil_contour_detection2(imgSmoothed, ...
%                                                   Px, Py, ...
%                                                   edgeThresh, ...
%                                                   round(radiiPupil), ...
%                                                   minFeatures);
% 
% %% Modified Original ellipse fitting method
% [~, inlierstest] = fit_ellipse_ransac(epxtest(:), epytest(:), radiiPupil + [-15, 15]);
% epx2test = epxtest(inlierstest);
% epy2test = epytest(inlierstest);
% ellipseResulttest = fit_ellipse(epx2test, epy2test);
% 
% pupiltest_x = ellipseResulttest.X0_in;
% pupiltest_y = ellipseResulttest.Y0_in;
% pupiltest_r1 = ellipseResulttest.a;
% pupiltest_r2 = ellipseResulttest.b;
% pupiltest_angle = -ellipseResulttest.phi;
% 
% ptestRadius = max([pupiltest_r1, pupiltest_r2]);
% ptestLx = pupiltest_r1*cos(phi)*cos(pupiltest_angle) - ...
%         sin(pupiltest_angle)*pupiltest_r2*sin(phi) + ...
%         pupiltest_x;
% ptestLy = pupiltest_r1*cos(phi)*sin(pupiltest_angle) + ...
%         cos(pupiltest_angle)*pupiltest_r2*sin(phi) + ...
%         pupiltest_y;
% ptestMx = pupiltest_x;
% ptestMy = pupiltest_y;



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