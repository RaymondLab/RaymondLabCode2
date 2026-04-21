function [trackParams] = setupEyeTracking_APP(app, trackParams)
%SETUPEYETRACKING_APP

% Define initial image to test
frame = 1;

% Detect pupil and CRs using set tracking parameters
ok = 0;
while ~ok && (frame < 10)

    % Mask image (with pupil contrast adjustment) based on pupil ROI
    [imgPupil] = readImg(trackParams.adjustmentValuesPupil, ...
                         trackParams.cam, ...
                         frame, ...
                         trackParams.maskPupil);

    % Mask image (with minimal contrasting) based on pupil ROI
    [img] = readImg(trackParams.defaultImgAdj, trackParams.cam, frame);

    % Set image geometry for pair-based CR tracker
    trackParams.crTrackParams.imgCenter = [size(img, 2)/2, size(img, 1)/2];
    trackParams.crTrackParams.isLeftCam = strcmp(trackParams.cameraSide, 'left');

    % Compute inter-CR distance (scale reference for spatial gating)
    trackParams.interCR_distance = sqrt( ...
        (trackParams.priMarker(1) - trackParams.secMarker(1))^2 + ...
        (trackParams.priMarker(2) - trackParams.secMarker(2))^2);

    % Initialize pair-based CR tracking state
    trackParams.pairState = struct( ...
        'prevPriPos',      trackParams.priMarker, ...
        'prevSecPos',      trackParams.secMarker, ...
        'expectedDist',    trackParams.interCR_distance, ...
        'afterBadGap',     false, ...
        'costSlidingWin',  zeros(trackParams.crTrackParams.slidingWindowSize, 1), ...
        'costWinCount',    0, ...
        'costWinIdx',      1, ...
        'distSlidingWin',  zeros(trackParams.crTrackParams.slidingWindowSize, 1), ...
        'distWinCount',    0, ...
        'distWinIdx',      1);

    % Detect pupil and CRs
    try
        [trackParams, frameData] = detectPupilCR_APP(app, img, imgPupil, [], [], [], trackParams, 1);
        ok = 1;
    catch msgid
        warning(msgid.message)
        ok = 0;
    end

    % Try next frame if bad image detection
    frame = frame + 1;
end

end