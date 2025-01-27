function [vars] = desaccadeMagnetVideo_APP(app, vars)
%DESACCADEMAGNETVIDEO_APP

%% Load data
loadAnalysisInfo_APP;

%% Prepare required data
freq = vars.stimFreq;
cutoffFreq = app.CutoffFrequencyEditField.Value;
windowPre = app.SaccadeWindow1EditField.Value;
windowPost = app.SaccadeWindow2EditField.Value;
minDataLength = app.MinimumGoodDataLengthEditField.Value;
mag1.saccadeThresh = app.Mag1ThresholdEditField.Value;
mag2.saccadeThresh = app.Mag2ThresholdEditField.Value;
vid.saccadeThresh = app.VidThresholdEditField.Value;


%% DESACCADE
if app.DesaccadeMethodDropDown.Value == "Squared Velocity Threshold"
    [mag1.saccades, ~, ~] = desaccadeSVT(mag1.pos_data_aligned, mag1.samplerate, 1, windowPre, windowPost, mag1.saccadeThresh, minDataLength, cutoffFreq);
    [mag2.saccades, ~, ~] = desaccadeSVT(mag2.pos_data_aligned, mag2.samplerate, 1, windowPre, windowPost, mag2.saccadeThresh, minDataLength, cutoffFreq);
    [vid.saccades, ~, ~]  = desaccadeSVT(vid.pos_data_upsampled_aligned, mag1.samplerate, 1, windowPre, windowPost, vid.saccadeThresh, minDataLength, cutoffFreq); 
elseif app.DesaccadeMethodDropDown.Value == "Moving Median MAD"
    [mag1.saccades, ~, ~, mag1.madThresh] = desaccadeMAD(mag1.pos_data_aligned, mag1.samplerate, 1, windowPre, windowPost, mag1.saccadeThresh, minDataLength, cutoffFreq);
    [mag2.saccades, ~, ~, mag2.madThresh] = desaccadeMAD(mag2.pos_data_aligned, mag2.samplerate, 1, windowPre, windowPost, mag2.saccadeThresh, minDataLength, cutoffFreq);
    [vid.saccades, ~, ~, vid.madThresh]  = desaccadeMAD(vid.pos_data_upsampled_aligned, mag1.samplerate, 1, windowPre, windowPost, vid.saccadeThresh, minDataLength, cutoffFreq); 
else
    error('"'+app.DesaccadeMethodDropDown.Value+'" is an invalid desaccading option.');
end

% Create the field "saccades_manual" if it does not already exist
if ~isfield(mag1, 'saccades_manual') || length(mag1.saccades_manual) ~= length(mag1.saccades)
    mag1.saccades_manual = false(length(mag1.saccades), 1);
end

if ~isfield(mag2, 'saccades_manual') || length(mag2.saccades_manual) ~= length(mag2.saccades)
    mag2.saccades_manual = false(length(mag2.saccades), 1);
end

if ~isfield(vid, 'saccades_manual') || length(vid.saccades_manual) ~= length(vid.saccades)
    vid.saccades_manual = false(length(vid.saccades), 1);
end

% Add video saccades to those found in magnets 1 and 2
mag1.saccades_all = mag1.saccades | mag1.saccades_manual | vid.saccades | vid.saccades_manual;
mag2.saccades_all = mag2.saccades | mag2.saccades_manual | vid.saccades | vid.saccades_manual;


%% SINE FIT AND SCALE FACTOR
% Magnet Channel 1 (saccades_all)
mag1Vel = mag1.vel_data_aligned;
mag1Vel(mag1.saccades_all) = nan;
[mag1.vel_amp, mag1.vel_phase, ~, mag1.vel_fit, mag1.vel_fitr2] = fit_sineWave(mag1Vel, mag1.samplerate, freq);

% Video (with mag1.saccades_all)
vidVel1 = vid.vel_data_upsampled_aligned;
vidVel1(mag1.saccades_all) = nan;
[vid.vel1_amp,  vid.vel1_phase,  ~, vid.vel1_fit, vid.vel1_fitr2] = fit_sineWave(vidVel1, mag1.samplerate, freq);

% Magnet Channel 1 Scaling Factor
scaleCh1 = abs(vid.vel1_amp / mag1.vel_amp);
mag1.vel_scale = mag1.align_sign*scaleCh1;

% Magnet Channel 2 (saccades_all)
mag2Vel = mag2.vel_data_aligned;
mag2Vel(mag2.saccades_all) = nan;
[mag2.vel_amp, mag2.vel_phase, ~, mag2.vel_fit, mag2.vel_fitr2] = fit_sineWave(mag2Vel, mag2.samplerate, freq);

% Video (with mag2.saccades_all)
vidVel2 = vid.vel_data_upsampled_aligned;
vidVel2(mag2.saccades_all) = nan;
[vid.vel2_amp,  vid.vel2_phase,  ~, vid.vel2_fit, vid.vel2_fitr2] = fit_sineWave(vidVel2, mag2.samplerate, freq);

% Magnet Channel 2
scaleCh2 = abs(vid.vel2_amp / mag2.vel_amp);
mag2.vel_scale = mag2.align_sign*scaleCh2;


%% Automatically select channel based on best r2 fit value
scaleCh1_Original = mag1.vel_scale;
scaleCh2_Original = mag2.vel_scale;

if mag1.vel_fitr2 > mag2.vel_fitr2
    scaleCh1 = scaleCh1_Original;
    scaleCh2 = 0.0;
    vid.vel_amp = vid.vel1_amp;
    vid.vel_fit = vid.vel1_fit;
    vid.vel_fitr2 = vid.vel1_fitr2;
    vid.saccades_all = mag1.saccades_all;
    vid.bestMagChannel = 1;
elseif mag2.vel_fitr2 > mag1.vel_fitr2
    scaleCh1 = 0.0;
    scaleCh2 = scaleCh2_Original;
    vid.vel_amp = vid.vel2_amp;
    vid.vel_fit = vid.vel2_fit;
    vid.vel_fitr2 = vid.vel2_fitr2;
    vid.saccades_all = mag2.saccades_all;
    vid.bestMagChannel = 2;
end

%% SAVE DATA
mag1Amp = mag1.vel_amp;
r2mag1 = mag1.vel_fitr2;
threshMagChan1 = mag1.saccadeThresh;
mag1Phase = mag1.vel_phase;

mag2Amp = mag2.vel_amp;
r2mag2 = mag2.vel_fitr2;
threshMagChan2 = mag2.saccadeThresh;
mag2Phase = mag2.vel_phase;

vidAmp = vid.vel_amp;
r2vid = vid.vel_fitr2;
threshVid = vid.saccadeThresh;

[~, filenameroot] = fileparts(cd);
save(fullfile(cd, [filenameroot,'.mat']), ...
     'scaleCh1', 'scaleCh2',...
     'vidAmp', 'mag1Amp', 'mag1Phase', 'mag2Amp', 'mag2Phase',...
     'r2mag1', 'r2mag2', 'r2vid', ...
     'threshVid', 'threshMagChan1', 'threshMagChan2', ...
     'freq');

saveAnalysisInfo_APP;

end