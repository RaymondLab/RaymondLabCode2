function [vars] = alignMagnetVideo_Linearity_APP(app, vars)
%ALIGNMAGNETVIDEO_LINEARITY_APP

%% Load data
% disp('Calculating Alignment Based on highest Linearity...')
loadAnalysisInfo_APP;

if isfield(vars.vid, 'bestMagChannel')
    mag1_scale = mag1.vel_scale;
    mag2_scale = mag2.vel_scale;
else
    mag1_scale = mag1.align_sign;
    mag2_scale = mag2.align_sign;
end

fc = app.CutoffFrequencyEditField.Value;
mag1_posfilt_data = butterworthfilter(mag1.pos_data, fc, mag1.samplerate);
mag1_posfilt_data = mag1_posfilt_data - nanmean(mag1_posfilt_data);
mag2_posfilt_data = butterworthfilter(mag2.pos_data, fc, mag1.samplerate);
mag2_posfilt_data = mag2_posfilt_data - nanmean(mag2_posfilt_data);
vid_posfilt_data_upsampled = butterworthfilter(vid.pos_data_upsampled, fc, mag1.samplerate);
vid_posfilt_data_upsampled = vid_posfilt_data_upsampled - nanmean(vid_posfilt_data_upsampled);

veltau = .01;
mag1.vel_data = movingslopeCausal(mag1_posfilt_data, round(mag1.samplerate*veltau))*mag1.samplerate;
mag2.vel_data = movingslopeCausal(mag2_posfilt_data, round(mag1.samplerate*veltau))*mag1.samplerate;
vid.vel_data_upsampled = movingslopeCausal(vid_posfilt_data_upsampled, round(mag1.samplerate*veltau)) * mag1.samplerate;


%% Calculate r2 values With Various Alignment Values
% Magnet Channel 1 Position
try
    [mag1.pos_linearity_r2, mag1.pos_linearity_maxr2, mag1.pos_linearity_maxr2Loc] = ...
        linearityAlignV2(mag1_posfilt_data*mag1_scale, vid_posfilt_data_upsampled);
catch
    disp('Error in magnet 1 position')
end

% Magnet Channel 1 Velocity
try
    [mag1.vel_linearity_r2, mag1.vel_linearity_maxr2, mag1.vel_linearity_maxr2Loc] = ...
        linearityAlignV2(mag1.vel_data*mag1_scale, vid.vel_data_upsampled);
catch
    disp('Error in magnet 1 velocity')

end

% Magnet Channel 2 Position
try
    [mag2.pos_linearity_r2, mag2.pos_linearity_maxr2, mag2.pos_linearity_maxr2Loc] = ...
        linearityAlignV2(mag2_posfilt_data*mag2_scale, vid_posfilt_data_upsampled); 
catch
    disp('Error in magnet 2 position')

end

% Magnet Channel 2 Velocity
try
    [mag2.vel_linearity_r2, mag2.vel_linearity_maxr2, mag2.vel_linearity_maxr2Loc] = ...
        linearityAlignV2(mag2.vel_data*mag2_scale, vid.vel_data_upsampled);
catch
    disp('Error in magnet 2 velocity')

end


%% Save data
saveAnalysisInfo_APP;
