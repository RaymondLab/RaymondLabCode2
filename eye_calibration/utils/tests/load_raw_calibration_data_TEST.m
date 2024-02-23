% Initialize Structs
vid = [];
mag1 = [];
mag2 = [];

% Load Camera Image Data
cam1FrameData_raw = load(fullfile(cd, 'videoresults_cam1.mat')).frameData;  % Camera 1
if size(cam1FrameData_raw,1) == 1
    cam1FrameData_raw = transpose(cam1FrameData_raw);
end

cam2FrameData_raw = load(fullfile(cd, 'videoresults_cam2.mat')).frameData;  % Camera 2
if size(cam2FrameData_raw,1) == 1
    cam2FrameData_raw = transpose(cam2FrameData_raw);
end

% Calculate Angular Head Position
vid.pos_raw = calceyeangle(cam1FrameData_raw, cam2FrameData_raw);
vid.percent_frames_dropped = sum(int64(isnan(vid.pos_raw)))*100 / length(vid.pos_raw);

% Import in Spike2 file
[~, filenameroot]= fileparts(cd);
fullfilename = fullfile(cd, [filenameroot,'.smr']);
smrfile = importSpike(fullfilename);