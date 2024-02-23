function [vars] = plotAlignedTraces_APP(app, vars)
%PLOTALIGNEDTRACES_APP

%% Align channels
vars = alignChannels_APP(vars, vars.chosenShiftVal);

%% Load data
loadAnalysisInfo_APP;

if isfield(vars.vid, 'bestMagChannel')
    mag1_scale = mag1.vel_scale;
    mag2_scale = mag2.vel_scale;
else
    mag1_scale = mag1.align_sign*vars.alignMagnetChannelScale;
    mag2_scale = mag2.align_sign*vars.alignMagnetChannelScale;
end

vid.pos_data_upsampled_aligned = cumtrapz(vid.vel_data_upsampled_aligned)./mag1.samplerate;
vid.pos_data_upsampled_aligned = vid.pos_data_upsampled_aligned - mean(vid.pos_data_upsampled_aligned);

mag1.vel_data_aligned_scaledInVel = mag1.vel_data_aligned.*mag1_scale;
mag2.vel_data_aligned_scaledInVel = mag2.vel_data_aligned.*mag2_scale;

vidTime = vid.time_upsampled_aligned;
mag1Time = mag1.time_aligned;
mag2Time = mag2.time_aligned;

mag1.pos_data_aligned_scaledInVel = cumtrapz(mag1.vel_data_aligned_scaledInVel)./mag1.samplerate;
mag2.pos_data_aligned_scaledInVel = cumtrapz(mag2.vel_data_aligned_scaledInVel)./mag2.samplerate;

% Center position at 2
mag1.pos_data_aligned_scaledInVel = mag1.pos_data_aligned_scaledInVel - mean(mag1.pos_data_aligned_scaledInVel);
mag2.pos_data_aligned_scaledInVel = mag2.pos_data_aligned_scaledInVel - mean(mag2.pos_data_aligned_scaledInVel);

if isfield(vars.vid, 'bestMagChannel')
    switch vars.vid.bestMagChannel
        case 1
            posYlim = 4*max([std(vid.pos_data_upsampled_aligned), std(mag1.pos_data_aligned_scaledInVel)]);
            velYlim = 4*max([std(vid.vel_data_upsampled_aligned), std(mag1.vel_data_aligned_scaledInVel)]);
        case 2
            posYlim = 4*max([std(vid.pos_data_upsampled_aligned), std(mag2.pos_data_aligned_scaledInVel)]);
            velYlim = 4*max([std(vid.vel_data_upsampled_aligned), std(mag2.vel_data_aligned_scaledInVel)]);
    end
else
    posYlim = 4*max([std(vid.pos_data_upsampled_aligned), std(mag1.pos_data_aligned_scaledInVel), std(mag2.pos_data_aligned_scaledInVel)]);
    velYlim = 4*max([std(vid.vel_data_upsampled_aligned), std(mag1.vel_data_aligned_scaledInVel), std(mag2.vel_data_aligned_scaledInVel)]);
end

%% Get cycle averages
[~, mag1_vinv_mean] = VOR_breakTrace(1000, 1, mag1.vel_data_aligned_scaledInVel);
[~, mag2_vinv_mean] = VOR_breakTrace(1000, 1, mag2.vel_data_aligned_scaledInVel);
[~, vid_mean]       = VOR_breakTrace(1000, 1, vid.vel_data_upsampled_aligned);

%% PLOT Position, Velocity, and Velocity Cycles
cla(app.UIAxesAlignedPositions);
cla(app.UIAxesAlignedVelocities);

% Position
plot(app.UIAxesAlignedPositions, vidTime, zeros(1,length(vidTime)), ':k', 'LineWidth',0.25, 'HandleVisibility','off');
hold(app.UIAxesAlignedPositions, 'on');
plot(app.UIAxesAlignedPositions, vidTime, vid.pos_data_upsampled_aligned, 'k', 'DisplayName','Video');
plot(app.UIAxesAlignedPositions, mag1Time, mag1.pos_data_aligned_scaledInVel, 'b', 'DisplayName','Magnet Channel 1');
plot(app.UIAxesAlignedPositions, mag2Time, mag2.pos_data_aligned_scaledInVel , 'r', 'DisplayName','Magnet Channel 2');
xlim(app.UIAxesAlignedPositions, [0, vidTime(end)]);
ylim(app.UIAxesAlignedPositions, [-posYlim, posYlim]);
ylabel(app.UIAxesAlignedPositions, 'Position (deg)');
legend(app.UIAxesAlignedPositions);
hold(app.UIAxesAlignedPositions, 'off');

% Velocity 
plot(app.UIAxesAlignedVelocities, vidTime, zeros(1,length(vidTime)), ':k', 'LineWidth',0.25, 'HandleVisibility','off');
hold(app.UIAxesAlignedVelocities, 'on');
plot(app.UIAxesAlignedVelocities, vidTime, vid.vel_data_upsampled_aligned, 'k', 'DisplayName','Video');
plot(app.UIAxesAlignedVelocities, mag1Time, mag1.vel_data_aligned_scaledInVel, 'b', 'DisplayName','Magnet Channel 1');
plot(app.UIAxesAlignedVelocities, mag2Time, mag2.vel_data_aligned_scaledInVel, 'r', 'DisplayName','Magnet Channel 2');
xlim(app.UIAxesAlignedVelocities, [0, vidTime(end)]);
ylim(app.UIAxesAlignedVelocities, [-velYlim, velYlim]);
xlabel(app.UIAxesAlignedVelocities, 'Time (s)');
ylabel(app.UIAxesAlignedVelocities, 'Velocity (deg/s)');
hold(app.UIAxesAlignedVelocities, 'off');

% Link the x-axes
linkaxes([app.UIAxesAlignedPositions, app.UIAxesAlignedVelocities], 'x');

% Velocity cycles
try
    cla(app.UIAxesAlignedVelocityCycles);
    velcYlim = 4*max([std(vid_mean), std(mag1_vinv_mean), std(mag2_vinv_mean)]);
    plot(app.UIAxesAlignedVelocityCycles, zeros(1,length(vid_mean)), ':k', 'LineWidth',0.25, 'HandleVisibility','off'); 
    hold(app.UIAxesAlignedVelocityCycles, 'on');
    plot(app.UIAxesAlignedVelocityCycles, vid_mean, 'k', 'DisplayName','Video');
    plot(app.UIAxesAlignedVelocityCycles, mag1_vinv_mean, 'b', 'DisplayName','Magnet Channel 1');
    plot(app.UIAxesAlignedVelocityCycles, mag2_vinv_mean, 'r', 'DisplayName','Magnet Channel 2');
    legend(app.UIAxesAlignedVelocityCycles);
    ylim(app.UIAxesAlignedVelocityCycles, [-velcYlim, velcYlim]);
    xlabel(app.UIAxesAlignedVelocityCycles, 'Time (ms)');
    ylabel(app.UIAxesAlignedVelocityCycles, 'Velocity (deg/s)');
    title(app.UIAxesAlignedVelocityCycles, 'Alignment of Cycle Averages');
    grid(app.UIAxesAlignedVelocityCycles, 'on');
catch
    % Do nothing
end


%% Save data
saveAnalysisInfo_APP;

end