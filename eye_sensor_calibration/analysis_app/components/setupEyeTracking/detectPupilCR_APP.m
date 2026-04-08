function [trackParams, newFrame, plotData] = detectPupilCR_APP(app, img, imgPupil, imgCRs, referenceFrame, newFrame, trackParams, plotOn)
% DETECTPUPILCR Detects pupil and corneal reflection in the eye image.
%
% Input:
% img        = input image with minimal contrasting
% imgPupil   = input image with pupil roi masking and contrasting
% imgCRs     = input image with CR roi masking and contrasting
% pupilStart = start point for starburst algorithm for pupil
% radiiPupil = guess of pupil radius
% edgeThresh = threshold for detecting pupil
%
% radiiCR    = guess of CR radius
% CRthresh   = threshold for CR
% CRfilter   = filter size for CR
%
% Output:
% pupil_ellipse = 5-vector of the ellipse parameters of pupil
%   [cx cy  a b theta]
%   cx - the x coordinate of ellipse center
%   cy - the y coordinate of ellipse center
%   a - the ellipse axis of x direction
%   b - the ellipse axis of y direction
%   theta - the orientation of ellipse
% cr1 and cr2 = 3-vector of the circle parameters of the corneal reflection
%   [crx cry crr]
%   crx - the x coordinate of circle center
%   cry - the y coordinate of circle center
%   crr - the radius of circle
% points        = actual points on pupil detected
% edgeThresh    = actual edgeThresh used for pupil
%
% Authors: Hannah Payne
% Date: 2014
%
%
%
% MODIFIED FROM:
%
% Starburst Algorithm
%
% This source code is part of the starburst algorithm.
% Starburst algorithm is free; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% Starburst algorithm is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with cvEyeTracker; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%
% Starburst Algorithm - Version 1.0.0
% Part of the openEyes ToolKit -- http://hcvl.hci.iastate.edu/openEyes
% Release Date:
% Authors : Dongheng Li <donghengli@gmail.com>
%           Derrick Parkhurst <derrick.parkhurst@hcvl.hci.iastate.edu>
% Copyright (c) 2005
% All Rights Reserved.

%% Process inputs
img = double(img);
imgPupil = double(imgPupil);

%% Preallocate
plotData.epx = nan;
plotData.epx2 = nan;
plotData.epy = nan;
plotData.epy2 = nan;
plotData.points = nan;

%% Find corneal reflections (v4 pipeline: top-hat + CHT + marker selection)

% Stage 1: Top-hat preprocessing on the base image
seTopHat = strel('disk', trackParams.topHatDiskRadius);
imgTopHat = double(imtophat(uint8(img), seTopHat));

% Stage 2: Circular Hough Transform on top-hat image
[~, circen, cirrad] = CircularHough_Grd(imgTopHat, trackParams.radiiCR, ...
    trackParams.CRthresh, trackParams.CRfilter, 1);

% Remove duplicate detections (30px threshold)
if size(circen, 1) > 1
    k = 1;
    while k < size(circen, 1)
        dups = sqrt(sum((repmat(circen(k,:), size(circen,1), 1) - circen).^2, 2)) < 30;
        dups(k) = false;
        circen(dups,:) = [];
        cirrad(dups) = [];
        k = k + 1;
    end
end

% Convert CHT circles to props struct for filtering/ranking
nCirc = size(circen, 1);
props = struct('WeightedCentroid',{}, 'Area',{}, 'Eccentricity',{}, ...
    'MaxIntensity',{}, 'Radius',{});
for ci = 1:nCirc
    cx = circen(ci, 1); cy = circen(ci, 2);
    cr = max(cirrad(ci), 1);
    props(ci).WeightedCentroid = [cx, cy];
    props(ci).Area = round(pi * cr^2);
    props(ci).Eccentricity = 0;
    iy = max(1, min(round(cy), size(img,1)));
    ix = max(1, min(round(cx), size(img,2)));
    props(ci).MaxIntensity = img(iy, ix);
    props(ci).Radius = cr;
end

% Stage 3: Candidate filtering by area and eccentricity
[validProps, ~] = filterCandidatesCR(props, ...
    trackParams.crMinArea, trackParams.crMaxArea, trackParams.crMaxEccentricity);

% Stage 4: Primary CR selection (soft-scored: prevPri + marker)
priFound = false; priWC = [NaN NaN]; priRad = NaN;
if ~isempty(validProps)
    bestPri = softScoreSelectCR(validProps, trackParams.priMarker, ...
        trackParams.prevPri, trackParams.wPrev, trackParams.wMarker);
    if ~isempty(bestPri)
        priWC = bestPri.WeightedCentroid;
        priRad = bestPri.Radius;
        priFound = true;
        trackParams.prevPri = priWC;
    end
end

% Find priIdx for excluding primary from secondary candidates
priIdx = 0;
if priFound
    for pidx = 1:numel(validProps)
        if isequal(validProps(pidx).WeightedCentroid, priWC)
            priIdx = pidx; break;
        end
    end
end

% Missing primary — use previous position
if ~priFound && ~isnan(trackParams.prevPri(1))
    priWC = trackParams.prevPri; priRad = NaN;
end

% Stage 5: Secondary CR selection (cone + soft-scored)
secWC = [NaN NaN]; secRad = NaN;
if ~isnan(priWC(1))
    if priIdx > 0
        secCandidates = validProps([1:priIdx-1, priIdx+1:end]);
    else
        secCandidates = validProps;
    end
    [secProp, ~, ~] = selectSecondaryCRfn( ...
        secCandidates, priWC, trackParams.cameraSide, ...
        trackParams.secConeHalfAngle, trackParams.secMarker, ...
        trackParams.prevSec, trackParams.wPrev, trackParams.wMarker);
    if ~isempty(secProp)
        secWC = secProp.WeightedCentroid;
        secRad = secProp.Radius;
        trackParams.prevSec = secWC;
    end
end
if all(isnan(secWC)) && ~isnan(trackParams.prevSec(1))
    secWC = trackParams.prevSec; secRad = NaN;
end

% Map primary/secondary to CR1(left)/CR2(right) for output compatibility
if strcmp(trackParams.cameraSide, 'left')
    crx = [priWC(1); secWC(1)];
    cry = [priWC(2); secWC(2)];
    crr_out = [priRad; secRad];
else
    crx = [secWC(1); priWC(1)];
    cry = [secWC(2); priWC(2)];
    crr_out = [secRad; priRad];
end

%% Remove corneal reflection from pupil image
if ~all(isnan(crx))
    totalMask = imgPupil > 204;  % Remove any bright spots
    InoCR = imgPupil;

    for i = 1:length(crx) % Remove each CR
        if isnan(crx(i)) || isnan(cry(i)), continue; end
        crMaskRadius = crr_out(i);
        if isnan(crMaskRadius), crMaskRadius = mean(trackParams.radiiCR); end
        [X, Y] = meshgrid(1:size(imgPupil,2), 1:size(imgPupil,1));
        maskCurr = ((X-crx(i)).^2 + (Y-cry(i)).^2) < crMaskRadius.^2;
        totalMask(maskCurr) = true;
    end

    totalMaskDilate = imdilate(totalMask,strel('disk', trackParams.spotMaskRadius));
    InoCR(totalMaskDilate) = NaN;

    gaussian_smooth_image = @(I, sigma) imfilter(I, fspecial('gaussian',[ceil(2.5*sigma),ceil(2.5*sigma)],sigma), 'symmetric');
    InoCR = gaussian_smooth_image(InoCR, 3);
else
    InoCR = imgPupil;
end

%% Find guess for pupil center using radial symmetry transform
if isempty(trackParams.pupilStart) || any(isnan(trackParams.pupilStart)) || trackParams.pupilStart(1)==0
    alphaFRST = 0.5;  % Sensitivity to radial symmetry
    imgRadialPupil = Radial_Sym_Transform(imgPupil, trackParams.radiiPupil, alphaFRST);
    [pupilY, pupilX] = find(min(imgRadialPupil(:))==imgRadialPupil);
    trackParams.pupilStart = [pupilX, pupilY];
end

%% Detect pupil borders using starburst algorithm
[epx, epy, ~] = starburst_pupil_contour_detection(InoCR, trackParams.pupilStart(1),...
    trackParams.pupilStart(2), trackParams.edgeThresh, round(trackParams.radiiPupil), trackParams.minfeatures);

[~, inliers] = fit_ellipse_ransac(epx(:), epy(:), trackParams.radiiPupil + [-15, 15]);

epx2 = epx(inliers);
epy2 = epy(inliers);

% Do better fit of resulting points
ellipseResult = fit_ellipse(epx2, epy2);

if isempty(ellipseResult.X0_in)
    newFrame.pupil_x = nan;
    newFrame.pupil_y = nan;
    newFrame.pupil_r1 = nan;
    newFrame.pupil_r2 = nan;
    newFrame.pupil_angle = nan;
else
    newFrame.pupil_x = ellipseResult.X0_in;
    newFrame.pupil_y = ellipseResult.Y0_in;
    newFrame.pupil_r1 = ellipseResult.a;
    newFrame.pupil_r2 = ellipseResult.b;
    newFrame.pupil_angle = -ellipseResult.phi;
end

newFrame.cr1_x = crx(1);
newFrame.cr1_y = cry(1);
newFrame.cr1_r = crr_out(1);
newFrame.cr2_x = crx(2);
newFrame.cr2_y = cry(2);
newFrame.cr2_r = crr_out(2);

points = [epx2(:), epy2(:)];

%% Add data to frameData object
plotData.epx = epx;
plotData.epx2 = epx2;
plotData.epy = epy;
plotData.epy2 = epy2;
plotData.points = points;

%% Plotting
if plotOn
   plotEyeTrackingImageFrame_APP(app, img, newFrame, plotData);
end
end

%% =====================================================================
%  LOCAL HELPER FUNCTIONS (v4 CR detection pipeline)
%  =====================================================================

function [validProps, rejectedProps] = filterCandidatesCR(props, minArea, maxArea, maxEcc)
% Filters CR candidates by area range and eccentricity.
    if isempty(props)
        validProps = props; rejectedProps = props; return;
    end
    keep = true(numel(props), 1);
    for k = 1:numel(props)
        if props(k).Area < minArea || props(k).Area > maxArea
            keep(k) = false; continue;
        end
        if props(k).Eccentricity > maxEcc
            keep(k) = false;
        end
    end
    validProps    = props(keep);
    rejectedProps = props(~keep);
end

function [secProp, inCone, outCone] = selectSecondaryCRfn( ...
    candidates, priWC, camSide, coneHalfAngle, secMarker, ...
    prevSec, wPrev, wMarker)
% Selects secondary CR using cone constraint from primary + soft scoring.
    secProp = []; inCone = []; outCone = [];
    if isempty(candidates) || any(isnan(priWC)), return; end

    if strcmp(camSide, 'left')
        coneDir = [1 0];
    else
        coneDir = [-1 0];
    end
    cosThresh = cosd(coneHalfAngle);
    isInCone  = false(numel(candidates), 1);

    for k = 1:numel(candidates)
        wc  = candidates(k).WeightedCentroid;
        vec = [wc(1) - priWC(1), wc(2) - priWC(2)];
        mag = sqrt(vec(1)^2 + vec(2)^2);
        if mag < eps, continue; end
        cosAngle = dot(vec, coneDir) / mag;
        if cosAngle >= cosThresh
            isInCone(k) = true;
        end
    end

    inCone  = candidates(isInCone);
    outCone = candidates(~isInCone);
    if isempty(inCone), return; end

    % Soft-scored selection among in-cone candidates
    secProp = softScoreSelectCR(inCone, secMarker, prevSec, wPrev, wMarker);
end

function bestProp = softScoreSelectCR(candidates, marker, prev, wPrev, wMarker)
% Selects the best candidate using a weighted score of distance to the
% previous position (temporal continuity) and the fixed marker (spatial
% prior). When no previous position exists, uses marker distance only.
    bestProp = [];
    if isempty(candidates), return; end

    scores = zeros(numel(candidates), 1);
    hasPrev = ~any(isnan(prev));

    for k = 1:numel(candidates)
        wc = candidates(k).WeightedCentroid;
        distToMarker = sqrt((wc(1)-marker(1))^2 + (wc(2)-marker(2))^2);
        if hasPrev
            distToPrev = sqrt((wc(1)-prev(1))^2 + (wc(2)-prev(2))^2);
            scores(k) = wPrev * distToPrev + wMarker * distToMarker;
        else
            scores(k) = distToMarker;
        end
    end

    [~, bestIdx] = min(scores);
    bestProp = candidates(bestIdx);
end