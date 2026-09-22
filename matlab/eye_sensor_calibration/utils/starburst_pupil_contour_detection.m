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

function [epx, epy, edge_thresh] = starburst_pupil_contour_detection(I, cx, cy, edge_thresh, radiiPupil, min_features)
%STARBURST_PUPIL_CONTOUR_DETECTION
% Input
% I = input image
% cx, cy = central start point of the feature detection process
% edge_thresh = best guess for the pupil contour threshold (30)

% Ouput
% epx = x coordinate of feature candidates [row vector]
% epy = y coordinate of feature candidates [row vector]

min_edge_thresh = 2;
N = 50;  % (Hannah = 100) Number of rays to use to detect feature points
min_candidate_features = N*min_features;  % Minimum number of pupil feature candidates
min_distance = radiiPupil(1)*.75;  % (Hannah = radiiPupil(1)*0.5) Distance from pupil center guess to start searching for edge
angle_spread = 100*pi/180;  % (Hannah = 60*pi/180)
min_change = 5;  % (Hannah = 6) Stop if sum of change in mean x and y over previous iteration is smaller than this (distance in pixels)
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
        % (Hannah) epx_init = epx;
        % (Hannah) epy_init = epy;
    end

    angle_normal = atan2(cy-epy, cx-epx);
    
    for i = 1:length(epx)
        [tepx, tepy, ~] = locate_edge_points(I, epx(i), epy(i), ...
            min_distance, ...  % (Hannah) min_distance*1.5
            angle_step*(edge_thresh/epd(i)), ...
            angle_normal(i), angle_spread, edge_thresh);
        epx = [epx tepx];
        epy = [epy tepy];
    end
    
    loop_count = loop_count+1;
    tcx(loop_count+1) = mean(epx);
    tcy(loop_count+1) = mean(epy);
    
    % (Hannah) sqrt((tcx(loop_count+1)-cx).^2 + (tcy(loop_count+1)-cy).^2) < min_change
    if abs(tcx(loop_count+1)-cx) + abs(tcy(loop_count+1)-cy) < min_change  % Threshold change in pupil position over loops % 10
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
epx1 = [];
epy1 = [];
dir1 = [];
ep_num = 0;  % ep stands for edge point
ep_num1 = 0;
step = 3;  % (Hannah = 2)
halfStep = step/2;
disMax = 40;
distances = dis:step:(dis+step*disMax);
angles = (angle_normal-angle_spread/2+0.0001):angle_step:(angle_normal+angle_spread/2);

distanceAmt = length(distances);
angleAmt = length(angles);

d_a_mat_x = round(cx + distances'*cos(angles));
d_a_mat_y = round(cy + distances'*sin(angles));

for a = 1:angleAmt
    
    if d_a_mat_y(1, a) > height || d_a_mat_y(1, a) < 1 || d_a_mat_x(1, a) > width || d_a_mat_x(1, a) < 1
        continue;
    end
    
    for d = 2:distanceAmt
        
        if d_a_mat_y(d,a) > height || d_a_mat_y(d,a) < 1 || d_a_mat_x(d,a) > width || d_a_mat_x(d,a) < 1
            break;
        end
        
        dw = I(d_a_mat_y(d,a),d_a_mat_x(d,a)) - I(d_a_mat_y(d-1,a),d_a_mat_x(d-1,a));
        
        if (dw >= edge_thresh)
            ep_num = ep_num+1;
            epx(ep_num) = d_a_mat_x(d-1,a)+halfStep;  % Edge point x coordinate
            epy(ep_num) = d_a_mat_y(d-1,a)+halfStep;  % Edge point y coordinate
            dir(ep_num) = dw;
            break;
        end 
    end
end
