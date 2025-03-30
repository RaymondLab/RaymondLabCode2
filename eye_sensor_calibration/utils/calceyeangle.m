function posH = calceyeangle(frameData_cam1, frameData_cam2, theta)

%% Angle between two cameras in degree  5/22 DONT CHANGE
if ~exist('theta','var')
    theta = 39.8;
end

%% Horizontal
d1 = [frameData_cam1.pupil_x] - [frameData_cam1.cr1_x];
d2 = [frameData_cam2.cr2_x] - [frameData_cam2.pupil_x];

posH = atand(sind(theta) ./ ((d1./d2) + cosd(theta)));
posH = posH - (theta/2);  % Re-centers posH to zero degrees