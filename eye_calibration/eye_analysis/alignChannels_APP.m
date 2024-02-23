function [vars] = alignChannels_APP(vars, lag)
%ALIGNCHANNELS_APP

%% Load data
loadAnalysisInfo_APP;


%% Add lag to channels
if lag < 0
    % Magnet channel 1
    mag1.pos_data_aligned = mag1.pos_data((-lag)+1:end);
    mag1.vel_data_aligned = mag1.vel_data((-lag)+1:end);

    % Magnet channel 2
    mag2.pos_data_aligned = mag2.pos_data((-lag)+1:end);
    mag2.vel_data_aligned = mag2.vel_data((-lag)+1:end);
    
    % Video
    vid.pos_data_upsampled_aligned = vid.pos_data_upsampled(1:end+lag);
    vid.vel_data_upsampled_aligned = vid.vel_data_upsampled(1:end+lag);
else

    % Magnet channel 1
    mag1.pos_data_aligned = mag1.pos_data(1:end-lag);
    mag1.vel_data_aligned = mag1.vel_data(1:end-lag);
    
    % Magnet channel 2
    mag2.pos_data_aligned = mag2.pos_data(1:end-lag);
    mag2.vel_data_aligned = mag2.vel_data(1:end-lag);

    % Video
    vid.pos_data_upsampled_aligned = vid.pos_data_upsampled((lag)+1:end);
    vid.vel_data_upsampled_aligned = vid.vel_data_upsampled((lag)+1:end);
end

%% If input vectors are not the same length, shorted the longer one
if length(mag1.pos_data_aligned) < length(vid.pos_data_upsampled_aligned)
    vid.pos_data_upsampled_aligned = vid.pos_data_upsampled_aligned(1:length(mag1.pos_data_aligned));
    vid.vel_data_upsampled_aligned = vid.vel_data_upsampled_aligned(1:length(mag1.pos_data_aligned));
elseif length(mag1.pos_data_aligned) > length(vid.pos_data_upsampled_aligned)
    % Magnet channel 1
    mag1.pos_data_aligned = mag1.pos_data_aligned(1:length(vid.pos_data_upsampled_aligned));
    mag1.vel_data_aligned = mag1.vel_data_aligned(1:length(vid.pos_data_upsampled_aligned));
    
    % Magnet channel 2
    mag2.pos_data_aligned = mag2.pos_data_aligned(1:length(vid.pos_data_upsampled_aligned));
    mag2.vel_data_aligned = mag2.vel_data_aligned(1:length(vid.pos_data_upsampled_aligned));
end

mag1.time_aligned = mag1.time(1:length(mag1.pos_data_aligned));
mag2.time_aligned = mag2.time(1:length(mag2.pos_data_aligned));
vid.time_upsampled_aligned = vid.time_upsampled(1:length(vid.pos_data_upsampled_aligned));

%% Save data
saveAnalysisInfo_APP;