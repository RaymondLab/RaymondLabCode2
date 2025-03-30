function vars = createAnalysisInfoData_APP(app, vars)

% Display a progress bar
pbar = uiprogressdlg(app.eyeCalibrationUIFigure, ...
    'Title','Please Wait', ...
    'Message','Loading and processing video and magnet channel data.');
pbar.Value = 0.01;

% Check for existing analysis info data
pathname = cd;
[~, filenameroot]= fileparts(pathname);
fullfilepath = fullfile(cd, [filenameroot '_analysis.mat']);
if isfile(fullfilepath)
    % Load data
    loadAnalysisInfo_APP;
    if isfield(vid, 'chosenShiftVal')
        vars.chosenShiftVal = vid.chosenShiftVal;
    end
else
    %% INITIAL SETUP
    try
        a = vars.mag1;
    catch
        vid = [];
        mag1 = [];
        mag2 = [];
    end
    freq = vars.stimFreq;
    
    %% LOAD VIDEOS
    % Camera 1
    pbar.Value = 0.1;
    try
        A = load(fullfile(cd, 'videoresults_cam1.mat'));
    catch
        error('Camera 1 Results Not Found: videoresults_cam1.mat')
    end
    rawFrameData_cam1 = A.frameData;
    
    % Camera 2
    pbar.Value = 0.2;
    try
        B = load(fullfile(cd, 'videoresults_cam2.mat'));
    catch
        error('Camera 2 Results Not Found: videoresults_cam2.mat')
    end
    rawFrameData_cam2 = B.frameData;
    
    % Convert pixel positions in video data to the angular eye position
    pbar.Value = 0.3;
    vid.pos_data = calceyeangle(rawFrameData_cam1, rawFrameData_cam2);
    vid.percent_frames_missed = sum(int64(isnan(vid.pos_data)))*100 / length(vid.pos_data);
    
    
    %% Sort out time stamps
    pbar.Value = 0.4;
    switch app.CameraTimeStampsDropDown.Value
        case "Camera 2"  % Camera 2 is default
            vid.time = [rawFrameData_cam1.time2];
        case "Camera 1"
            vid.time = [rawFrameData_cam1.time1];
        case "Median"
            % Median Times between frames
            vid.time = median([rawFrameData_cam1.time1; rawFrameData_cam1.time2]);
    end
    vid.time = vid.time - vid.time(1);
    
    pbar.Value = 0.5;
    [tscale, ~] = fminsearchbnd(@(x) vidTimeFcn_APP(app,vid.time,vid.pos_data',freq,x), 1, .7, 1.4);  % old
    vid.time = vid.time / tscale;  % .995; % old
    vid.samplerate = 1 / mean(diff(vid.time));  % old
    
    
    %% SETUP MAGNET CHANNEL
    % Load Magnet Data
    pbar.Value = 0.6;
    [~, filenameroot]= fileparts(cd);
    fullfilename = fullfile(cd, [filenameroot,'.smr']);
    
    try
        if app.EyeMagnetLocationDropDown.Value == "Left"
            % Left Eye (default)
            rawMagnetData = importSpike(fullfilename, [4 5 6 10 31]);
        else
            % Right Eye
            rawMagnetData = importSpike(fullfilename, [4 7 8 10 31]);
        end
        
        % Select Proper Magnet window/segment
        lightpulses = rawMagnetData(end).data;
        segmentStart = lightpulses(1);
        segmentEnd = segmentStart+vid.time(end);
        
        rawMagnetData = resettime(datseg(rawMagnetData, [segmentStart segmentEnd]));
        
        mag1.pos_data = double(rawMagnetData(2).data);
        mag1.samplerate = rawMagnetData(2).samplerate;
        mag1.time = dattime(rawMagnetData(2));
        
        mag2.pos_data = double(rawMagnetData(3).data);
        mag2.samplerate = rawMagnetData(3).samplerate;
        mag2.time = dattime(rawMagnetData(3));
    catch
        rawMagnetData = importSpike(fullfilename, [2 7 8 11]);
        
        % Select Proper Magnet window/segment
        lightpulses = rawMagnetData(end).data;
        segmentStart = lightpulses(1);
        segmentEnd = segmentStart+vid.time(end);
        
        rawMagnetData = resettime(datseg(rawMagnetData,[segmentStart segmentEnd]));
        
        mag1.pos_data = double(rawMagnetData(2).data);
        mag1.samplerate = rawMagnetData(2).samplerate;
        mag1.time = dattime(rawMagnetData(2));
        
        mag2.pos_data = double(rawMagnetData(3).data);
        mag2.samplerate = rawMagnetData(3).samplerate;
        mag2.time = dattime(rawMagnetData(3));
    end
    
    % Alignment (positive/negative) sign for magnet channels
    mag1.align_sign = 1.0;
    mag2.align_sign = 1.0;
    
    
    %% Upsample Video Traces (OLD METHOD)
    pbar.Value = 0.8;
    vid.pos_data_upsampled = interp1(vid.time, vid.pos_data, mag1.time(:), 'linear');
    vid.pos_data_upsampled = inpaint_nans(vid.pos_data_upsampled);
    vid.time_upsampled = mag1.time;

end

%% SAVE DATA
pbar.Value = 0.9;
saveAnalysisInfo_APP;

close(pbar);
pause(0.01);

end