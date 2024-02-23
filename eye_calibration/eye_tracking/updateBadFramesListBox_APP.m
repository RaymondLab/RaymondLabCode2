function [vars, frameData] = updateBadFramesListBox_APP(app, vars, trackParams, frameData, nanBadFrames)
%UPDATEBADFRAMESLISTBOX_APP Updates bad frames list box on app.

if ~exist('nanBadFrames', 'var')
    nanBadFrames = 0;
end

cam_badFrames = findBadFrames(frameData);
badFrames = num2cell(find(cam_badFrames));
badFrames = cellfun(@num2str, badFrames, 'UniformOutput',false);

if (length(badFrames) == 1) && (str2double(badFrames{1}) == trackParams.nImages)
    badFrames = {'None'};
    app.RedoSelectedFrameButton.Enable = "off";
    app.RedoCustomFrameButton.Enable = "on";
else
    badFrames(ismember(badFrames,num2str(trackParams.nImages))) = [];
    app.RedoSelectedFrameButton.Enable = "on";
    app.RedoCustomFrameButton.Enable = "on";
end

if trackParams.cam == 1
    vars.cam1_badFrames = cam_badFrames;
    vars.cameraHasTrackingData(1) = 1;
else
    vars.cam2_badFrames = cam_badFrames;
    vars.cameraHasTrackingData(2) = 1;
end

app.BadFramesListBox.Items = badFrames;

if nanBadFrames && ~isnan(badFrames{1})
    for i = 1:length(badFrames)
        badFrame = str2double(badFrames{i});
        frameData(badFrame).cr1_x = nan;
        frameData(badFrame).cr1_y = nan;
        frameData(badFrame).cr1_r = nan;

        frameData(badFrame).cr2_x = nan;
        frameData(badFrame).cr2_y = nan;
        frameData(badFrame).cr2_r = nan;

        frameData(badFrame).pupil_x = nan;
        frameData(badFrame).pupil_y = nan;
        frameData(badFrame).pupil_r1 = nan;
        frameData(badFrame).pupil_r2 = nan;
        frameData(badFrame).pupil_angle = nan;
    end
end

end