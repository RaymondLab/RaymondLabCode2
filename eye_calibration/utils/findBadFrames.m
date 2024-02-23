function badFrames = findBadFrames(frameData)
%FINDBADFRAMES Finds the indices of frames with poor detections.

threshCR = 30;
threshpupilr = 10;

%% Find NAN frames
frameData_table = struct2table(frameData);
frameData_nanLocations = ismissing(frameData_table);
badFrames = any(frameData_nanLocations, 2);

%% Find bad corneal reflection detections
cr1_x_deriv = diff([frameData.cr1_x]);
cr1_y_deriv = diff([frameData.cr1_y]);
cr2_x_deriv = diff([frameData.cr2_x]);
cr2_y_deriv = diff([frameData.cr2_y]);

cr1_distanceChange = sqrt(cr1_x_deriv.^2 + cr1_y_deriv.^2);
cr2_distanceChange = sqrt(cr2_x_deriv.^2 + cr2_y_deriv.^2);

bad_cr1 = cr1_distanceChange > threshCR;
bad_cr2 = cr2_distanceChange > threshCR;

bad_cr1 = [bad_cr1, 0];
bad_cr2 = [bad_cr2, 0];

%% Find bad pupil detections
pupil_r1_deriv = abs(diff([frameData.pupil_r1]));
pupil_r2_deriv = abs(diff([frameData.pupil_r2]));

bad_pupil1 = pupil_r1_deriv > threshpupilr;
bad_pupil2 = pupil_r2_deriv > threshpupilr;

bad_pupil1 = [bad_pupil1, 0];
bad_pupil2 = [bad_pupil2, 0];

%% Add bad CRs list of bad frames
badFrames = badFrames' | bad_cr1 | bad_cr2 | bad_pupil1 | bad_pupil2;






