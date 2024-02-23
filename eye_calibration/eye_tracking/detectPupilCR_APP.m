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

%% Process imputs
img = double(img);
imgCRs = double(imgCRs);
imgPupil = double(imgPupil);

%% Preallocate
plotData.epx = nan;
plotData.epx2 = nan;
plotData.epy = nan;
plotData.epy2 = nan;
plotData.points = nan;

%% Find corneal reflections
% only search through portion of image, Corneal reflections shouldn't move
try 
    if ~isempty([referenceFrame.cr1_x,referenceFrame.cr2_x]) && ~any(isnan([referenceFrame.cr1_x,referenceFrame.cr2_x]))

        TOP    = floor(max(min([referenceFrame.cr1_y,referenceFrame.cr2_y]) - ...
                           max([referenceFrame.cr1_r,referenceFrame.cr2_r])*2, 1 ));
        BOTTOM = floor(min(max([referenceFrame.cr1_y,referenceFrame.cr2_y]) + ...
                           max([referenceFrame.cr1_r,referenceFrame.cr2_r])*2, size(imgCRs,1) ));
        LEFT   = floor(max(min([referenceFrame.cr1_x,referenceFrame.cr2_x]) - ...
                           max([referenceFrame.cr1_r,referenceFrame.cr2_r])*2, 1 ));
        RIGHT  = floor(min(max([referenceFrame.cr1_x,referenceFrame.cr2_x]) + ...
                           max([referenceFrame.cr1_r,referenceFrame.cr2_r])*2, size(imgCRs,2) ));

        newImg = imgCRs(TOP:BOTTOM, LEFT:RIGHT);
        [~, circen, crr] = CircularHough_Grd(newImg, trackParams.radiiCR, trackParams.CRthresh, trackParams.CRfilter, 1);

        circen(:,1) = circen(:,1) + LEFT-1;
        circen(:,2) = circen(:,2) + TOP-1;
        
        for i = length(circen):-1:1
            if imgCRs(floor(circen(i,2)), floor(circen(i,1))) == 0
                circen(i,:) = [];
                crr(i) = [];
            elseif imgCRs(floor(circen(i,2)), floor(circen(i,1))) < 150
                circen(i,:) = [];
                crr(i) = [];
            end
        end
    else
        [A, circen, crr] = CircularHough_Grd(imgCRs, trackParams.radiiCR, trackParams.CRthresh, trackParams.CRfilter, 1);
    end
catch
    [A, circen, crr] = CircularHough_Grd(imgCRs, trackParams.radiiCR, trackParams.CRthresh, trackParams.CRfilter, 1);
end

maxvals = diag(imgCRs(round(circen(:,2)),round(circen(:,1))));

%% Check for duplicates
i=1;
while i<size(circen,1)
    duplicates = sqrt(sum((repmat(circen(i,:),size(circen,1),1) - circen).^2,2)) <30;
    duplicates(i)=0;
    circen(duplicates,:)=[];
    crr(duplicates) = [];
    maxvals(duplicates) = [];
    i=i+1;
end

if size(circen,1)>4
    circen = circen(1:4,:);
end

cry = circen(:,2);
crx = circen(:,1);
[cry, I] = sort(cry,1,'descend');
crr = crr(I);
crx = crx(I);

%% Make sure first 2 crs go from left to right
if size(circen,1)>=2
    [crx(1:2), I] = sort(crx(1:2));
    cry(1:2) = cry(I);
    crr(1:2) = crr(I);
else
    cry(2) = NaN;
    crx(2) = NaN;
    crr(2) = NaN;
end

%% Remove corneal reflection
if ~isempty(circen)
    totalMask = imgPupil > 204;  % Remove any bright spots
    InoCR = imgPupil;
    
    for i = 1:length(crx) % Remove each CR
        [X, Y] = meshgrid(1:size(imgPupil,2), 1:size(imgPupil,1));
        maskCurr = ((X-crx(i)).^2 + (Y-cry(i)).^2) < crr(i).^2;
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
newFrame.cr1_r = crr(1);
newFrame.cr2_x = crx(2);
newFrame.cr2_y = cry(2);
newFrame.cr2_r = crr(2);

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