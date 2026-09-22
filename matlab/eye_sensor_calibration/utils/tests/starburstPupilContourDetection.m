% Starburst Algorithm 
%
% This function implements the Starburst Algorithm, a feature-based approach
% for detecting the contour of the pupil in eye-tracking applications. It uses
% a radial search pattern emanating from an initial guess point to locate edge
% points that are likely part of the pupil boundary. The algorithm iteratively
% refines these points to converge on the pupil contour.
%
% Usage:
% [epx, epy, edge_thresh, epx_init, epy_init] = starburstPupilContourDetection(I, cx, cy, edge_thresh, radii_pupil, min_features, min_edge_thresh)
%
% Inputs:
% - I: Input image (grayscale) in which pupil contour is to be detected.
% - cx, cy: Initial guess for the pupil center (x and y coordinates).
% - edge_thresh: Initial threshold for edge detection. The algorithm adjusts this
%                threshold dynamically based on feature detection. Default: 30
% - radii_pupil: Expected radii range of the pupil. Only the first value is used. It is scaled 
%               to set the distance at which to start searching for an edge.
% - min_features: Minimum fraction of total features (relative to the number of rays)
%                required to proceed with the detection.
% - min_edge_thresh: Minimum allowable edge threshold to limit the adjustment of
%                    edge_thresh during detection. Default: 2
%
% Outputs:
% - epx, epy: Coordinates (x and y) of the detected feature points that likely
%             belong to the pupil contour.
% - edge_thresh: The final edge threshold value used for detecting the feature points.
% - epx_init, epy_init: Initial sets of feature points detected before iterative refinement.
%
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
%           Edited by Hannah Payne, 2021
% Copyright (c) 2005
% All Rights Reserved.

function [epx, epy, edge_thresh, epx_init, epy_init] = starburstPupilContourDetection(I, cx, cy, edge_thresh, radii_pupil, min_features, min_edge_thresh)

if ~exist('min_edge_thresh', 'var') || isempty(min_edge_thresh)
    min_edge_thresh = 2;
end

if isempty(edge_thresh) || edge_thresh<= min_edge_thresh
    edge_thresh = 30;  % edge_threshold = best guess for the pupil contour threshold (30)
end

N = 100;  % Number of rays to use to detect feature points. Edited from 50
min_candidate_features = N*min_features;  % Minimum number of pupil feature candidates
min_distance = radii_pupil(1)*.5;  % Distance from pupil center guess to start searching for edge 
angle_spread = 60*pi/180;               
min_change = 6;  % Stop if sum of change in mean x and y over previous iteration is smaller than this (distance in pixels)
loop_count = 0;
max_loop_count = 10;
tcx(loop_count+1) = cx;
tcy(loop_count+1) = cy;
angle_step= 2*pi/N;

while edge_thresh > min_edge_thresh && loop_count <= max_loop_count
    epx = [];
    epy = [];
    while length(epx) < min_candidate_features && edge_thresh > min_edge_thresh
        [epx, epy, epd] = locate_edge_points(I, cx, cy, min_distance, angle_step, 0, 2*pi, edge_thresh);
        if length(epx) < min_candidate_features
            edge_thresh = edge_thresh - 1;
        end
        epx_init = epx;
        epy_init = epy;
    end

    angle_normal = atan2(cy-epy, cx-epx);
    
    for i = 1:length(epx)
        [tepx, tepy, ~] = locate_edge_points(I, epx(i), epy(i), ...
            min_distance*1.5, angle_step*(edge_thresh/epd(i)), ...
            angle_normal(i), angle_spread, edge_thresh);
        epx = [epx tepx];
        epy = [epy tepy];
    end
    
    loop_count = loop_count+1;
    tcx(loop_count+1) = mean(epx);
    tcy(loop_count+1) = mean(epy);
    if sqrt((tcx(loop_count+1)-cx).^2 + (tcy(loop_count+1)-cy).^2) < min_change % threshold change in pupil position over loops
        break;
    end
    cx = mean(epx);
    cy = mean(epy);
end


function [epx, epy, dir] = locate_edge_points(I, cx, cy, dis, angle_step, angle_normal, angle_spread, edge_thresh)

[height, width] = size(I);
epx = [];
epy = [];
dir = [];
ep_num = 0;  % ep stands for edge point
step = 2; 
maxDist = 2;
maxstep = round(dis*maxDist/step);
dists = dis + step*(0:maxstep)';

for angle=(angle_normal-angle_spread/2+0.0001):angle_step:(angle_normal+angle_spread/2)
    cosangle = cos(angle);
    sinangle = sin(angle);
    
    % Generate candidate points
    ps = [round(cx + dists*cosangle)  round(cy + dists*sinangle)]; 
    
    % Make sure points are within image bounds
    mask =  ps(:,2) < height & ps(:,2) > 1 & ps(:,1) < width & ps(:,1) > 1;
    ps = ps(mask,:);
        
    inds = ps(:,2) + (ps(:,1) - 1).*height; 
    vals = I(inds);    
    
    d = diff(vals);    
    ind = find(d>edge_thresh,1);
    
    if ~isempty(ind)
        ep_num = ep_num+1;
        epx(ep_num) = ps(ind,1)+step/2;  % Edge point x coordinate
        epy(ep_num) = ps(ind,2)+step/2;  % Edge point y coordinate
        dir(ep_num) = d(ind);        
    end
end