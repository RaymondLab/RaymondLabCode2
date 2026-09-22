function [vars, trackParams, frameData] = loadExistingTrackingData_APP(app, vars, trackParams, frameData)
%RESTOREPREVIOUSFRAMEDATA_APP Restores previous eye tracking analysis data.

if ~vars.isValidFolder || ~vars.cameraHasTrackingData(trackParams.cam)
    frameData = [];
    app.BadFramesListBox.Items = {'None'};
    app.RedoSelectedFrameButton.Enable = "off";
    app.RedoCustomFrameButton.Enable = "off";
    app.RunEyeTrackingButton.Enable = "off";
    cla(app.UIAxesTab1EyeTracking);
    cla(app.UIAxesTab1PupilPos);
    cla(app.UIAxesTab1CRPos);
    cla(app.UIAxesTab1Radii);
    return
end

if ~exist('frameData', 'var') || isempty(frameData)
    % Load previous eye tracking analysis data
    load(['videoresults_cam',num2str(trackParams.cam),'.mat'], 'frameData');
end

% Ensure frameData has the correct dimensions, i.e. (nFrames,1)
% This is a fix for D241B data which seems to all be reversed
if size(frameData,1) == 1
    frameData = transpose(frameData);
    save(['videoresults_cam',num2str(trackParams.cam),'.mat'], 'frameData');
end

% Load last frame from images and apply processing
trackParams.nImages = max(size(frameData));
lastIdx = round(trackParams.nImages - 1);
img = readImg(trackParams.defaultImgAdj, trackParams.cam, lastIdx);

% Plot Eye Tracking
plotEyeTrackingImageFrame_APP(app, img, frameData(lastIdx));

% Plot pupil,CR1, and CR2 positions and radii
plotEyeTrackingPositionsRadii_APP(app, frameData);

% Plot a scrollbar and a listener
sb = drawline(app.UIAxesTab1PupilPos, ...
              'Position',[lastIdx -100;lastIdx 100], ...
              'Color','g', ...
              'Label',sprintf(['Frame ', num2str(lastIdx)]), ...
              'InteractionsAllowed','translate');
addlistener(sb, 'MovingROI', @(src,evt) scrollbarMoved_APP(src,evt,app,trackParams,frameData));

% Update bad frames list on app
[vars, frameData] = updateBadFramesListBox_APP(app, vars, trackParams, frameData);

app.RunEyeTrackingButton.Enable = "on";

end