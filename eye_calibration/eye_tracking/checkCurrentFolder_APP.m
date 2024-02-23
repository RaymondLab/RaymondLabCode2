function [vars] = checkCurrentFolder_APP(app, vars)
%CHECKCURRENTFOLDER_APP Checks Current Folder for required/existing data.

% Check whether Current Folder has required files for calibration
if ~isfile('img1.tiff') || ~isfile('img2.tiff') || ~isfile('time.mat')
    vars.isValidFolder = 0;
    uialert(app.eyeCalibrationUIFigure, ...
        'Current Folder is not a valid calibration folder.', ...
        'Invalid Folder');
else
    vars.isValidFolder = 1;
end

% Check whether folder already has existing eye tracking data
vars.cameraHasTrackingData = [isfile('videoresults_cam1.mat'), ...
                              isfile('videoresults_cam2.mat')];

end