function [vars] = desaccadePreprocessing_APP(app, vars)
%DESACCADEPREPROCESSING_APP

% Load data
loadAnalysisInfo_APP;

% Preprocess required data
mag1_posfilt_data_aligned = butterworthfilter(mag1.pos_data_aligned, app.CutoffFrequencyEditField.Value, mag1.samplerate);
mag2_posfilt_data_aligned = butterworthfilter(mag2.pos_data_aligned, app.CutoffFrequencyEditField.Value, mag1.samplerate);
vid_posfilt_data_aligned  = butterworthfilter(vid.pos_data_upsampled_aligned, app.CutoffFrequencyEditField.Value, mag1.samplerate);

veltau = .01;
mag1.vel_data_aligned          = movingslopeCausal(mag1_posfilt_data_aligned, round(mag1.samplerate*veltau)) * mag1.samplerate;
mag2.vel_data_aligned          = movingslopeCausal(mag2_posfilt_data_aligned, round(mag1.samplerate*veltau)) * mag1.samplerate;
vid.vel_data_upsampled_aligned = movingslopeCausal(vid_posfilt_data_aligned, round(mag1.samplerate*veltau)) * mag1.samplerate;

% Save data
saveAnalysisInfo_APP;

end