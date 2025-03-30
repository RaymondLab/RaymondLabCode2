function [trackParams] = setupEyeTracking_APP(app, trackParams)
%SETUPEYETRACKING_APP

% Define initial image to test
frame = 1;

% Detect pupil and CRs using set tracking parameters
ok = 0;
while ~ok && (frame < 10)

    % Mask image (with CR contrast adjustment) based on CR ROI
    [imgCRs] = readImg(trackParams.adjustmentValuesCRs, ...
                       trackParams.cam, ...
                       frame, ...
                       trackParams.maskCRs);

    % Mask image (with pupil contrast adjustment) based on pupil ROI
    [imgPupil] = readImg(trackParams.adjustmentValuesPupil, ...
                         trackParams.cam, ...
                         frame, ...
                         trackParams.maskPupil);

    % Mask image (with minimal contrasting) based on pupil ROI
    [img] = readImg(trackParams.defaultImgAdj, trackParams.cam, frame);

    % Detect pupil and CRs
    try
        [trackParams, frameData] = detectPupilCR_APP(app, img, imgPupil, imgCRs, [], [], trackParams, 1);
        ok = 1;
    catch msgid
        warning(msgid.message)
        ok = 0;
    end

    % Try next frame if bad image detection
    frame = frame + 1;
end

end