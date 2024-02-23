function [trackParams, frameData] = eyeAnalysis_APP(app, trackParams)

%% Parameters & Prep
theta = 0:.1:2*pi;
phi = linspace(0, 2*pi, 100);
lastGoodFrame = 1;

%% Get Images
% Display a progress bar
pbar = uiprogressdlg(app.eyeCalibrationUIFigure, ...
    'Title','Please Wait', ...
    'Message','Loading images from .tiff file.');
pbar.Value = 0.1;
[imgStack, ~] = getImageStack(['img', num2str(trackParams.cam), '.tiff']);
trackParams.nImages = size(imgStack, 3);

%% Pre-process Images
pbar.Value = .4; 
pbar.Message = 'Preprocessing images.';
imgStackCrs = preprocessImages(imgStack, trackParams.adjustmentValuesCRs, trackParams.maskCRs);
pbar.Value = .6; 
imgStackPupil = preprocessImages(imgStack, trackParams.adjustmentValuesPupil, trackParams.maskPupil);
pbar.Value = .8; 
imgStack = preprocessImages(imgStack, trackParams.defaultImgAdj);

pbar.Value = .9; 
%% Load results time file & Preallocate frameData Object
load(fullfile(cd,'time.mat'));

n = length(results.time1);
frameData = struct('cr1_x', nan, 'cr1_y', nan, 'cr1_r', nan, ...
                   'cr2_x', nan, 'cr2_y', nan, 'cr2_r', nan, ...
                   'pupil_x', nan, 'pupil_y', nan, ...
                   'pupil_r1', nan, 'pupil_r2', nan, ...
                   'pupil_angle', nan, ...
                   'time1', num2cell(results.time1), ...
                   'time2', num2cell(results.time2), ...
                   'time3', num2cell(results.time3));

pbar.Value = 0.97;
pause(0.01);
close(pbar);
               
%% Start looping
for i = 1:n-1
    
    %% every 15 frames, update position and radii figures
    if mod(i,15) == 0
        if i == 15
            hp = plotEyeTrackingPositionsRadii_APP(app, frameData);
        else
            hp = plotEyeTrackingPositionsRadii_APP(app, frameData, hp);
        end
    end
    
    %% Load image
    img = imgStack(:,:,i);
    imgCRs = imgStackCrs(:,:,i);
    imgPupil = imgStackPupil(:,:,i);
    
    %% Use most recent good frame as a starting point
    if i ~= 1
        lastGoodFrame = i;
        
        while (any(structfun(@isnan, frameData(lastGoodFrame))) || any(structfun(@isempty, frameData(lastGoodFrame)))) && lastGoodFrame > 0
            lastGoodFrame = lastGoodFrame - 1;
        end

        maxRadii = nanmax(frameData(lastGoodFrame).pupil_r1, frameData(lastGoodFrame).pupil_r2);
        trackParams.radiiPupil(2) = round(maxRadii*1.15);
    end
    
    %% Detect Pupil  
    try
        [trackParams, frameData(i), plotData] = detectPupilCR_APP(app, img, imgPupil, imgCRs, frameData(lastGoodFrame), frameData(i), trackParams, 0);
    catch msgid
        fprintf('Error in img %i:\n', i)
        disp(msgid.message)
        trackParams.edgeThresh = 35;
    end
    
    %% Plotting
    if i == 1  %CHECK IF WE CAN ADD PLOTDATA TO SAVE FILE?
        plots = plotEyeTrackingImageFrame_APP(app, img, frameData(i), plotData);
    else
        plots = plotEyeTrackingImageFrame_APP(app, img, frameData(i), plotData, plots);
    end
    pause(0.01)  % Needed for things to actually plot.
end

%% Plot a scrollbar and a listener
lastIdx = round(trackParams.nImages - 1);
sb = drawline(app.UIAxesTab1PupilPos, ...
              'Position',[lastIdx -100;lastIdx 100], ...
              'Color','g', ...
              'Label',sprintf(['Frame ', num2str(lastIdx)]), ...
              'LabelVisible','hover', ...
              'InteractionsAllowed','translate');
addlistener(sb, 'MovingROI', @(src,evt) scrollbarMoved_APP(src,evt,app,trackParams,frameData));

%% Save Data
% Ensure frameData has the correct dimensions, i.e. (nFrames,1)
% This is a fix for D241B data which seems to all be reversed
if size(frameData,1) == 1
    frameData = transpose(frameData);
end
save(['videoresults_cam', num2str(trackParams.cam), '.mat'], 'frameData');
fprintf('\nRun Eye Tracking results saved!');

%% Save Plot Summary Pdf and Fig
saveplotEyeTrackingCam(frameData, trackParams.cam);