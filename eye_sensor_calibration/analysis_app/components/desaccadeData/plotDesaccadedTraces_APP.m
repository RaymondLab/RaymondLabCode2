function [vars] = plotDesaccadedTraces_APP(app, vars, updateOnly)
%PLOTDESACCADEDTRACES_MAD_APP

%% Load data
loadAnalysisInfo_APP;

%% Prepare data for plotting
if app.DesaccadeMethodDropDown.Value == "Squared Velocity Threshold"
    % Squared Velocity Error Relative to Sinusoidal Fit
    mag1VelThreshSac = (mag1.vel_data_aligned - mag1.vel_fit').^2;
    mag2VelThreshSac = (mag2.vel_data_aligned - mag2.vel_fit').^2;
    vidVelThreshSac  = (vid.vel_data_upsampled_aligned - vid.vel_fit').^2;
    % Threshold Arrays for Plotting
    mag1Thresh = ones(length(mag1VelThreshSac),1).*mag1.saccadeThresh;
    mag2Thresh = ones(length(mag2VelThreshSac),1).*mag2.saccadeThresh;
    vidThresh = ones(length(vidVelThreshSac),1).*vid.saccadeThresh;
    % Calculate the Y-Axis Limits for the Threshold Plot
    mag1ThreshYLim = mag1.saccadeThresh + (mag1.saccadeThresh*0.5);
    mag2ThreshYLim = mag2.saccadeThresh + (mag2.saccadeThresh*0.5);
    vidThreshYLim  = vid.saccadeThresh + (vid.saccadeThresh*0.5);
    threshYLabel = 'Squared Velocity (V2/s2)';
elseif app.DesaccadeMethodDropDown.Value == "Moving Median MAD"
    % Velocity Error Relative to Sinusoidal Fit
    mag1MovMedVel = (mag1.vel_data_aligned - mag1.vel_fit');
    mag2MovMedVel = (mag2.vel_data_aligned - mag2.vel_fit');
    vidMovMedVel  = (vid.vel_data_upsampled_aligned - vid.vel_fit');
    % Absolute Moving Zero-Median of the Velocity Error
    mag1VelThreshSac = abs(mag1MovMedVel - movmedian(mag1MovMedVel, 1000));
    mag2VelThreshSac = abs(mag2MovMedVel - movmedian(mag2MovMedVel, 1000));
    vidVelThreshSac = abs(vidMovMedVel - movmedian(vidMovMedVel, 1000));
    % Threshold Arrays for Plotting
    mag1Thresh = mag1.madThresh;
    mag2Thresh = mag2.madThresh;
    vidThresh = vid.madThresh;
    % Calculate the Y-Axis Limits for the Threshold Plot
    mag1ThreshYLim = 3*std(mag1.madThresh);
    mag2ThreshYLim = 3*std(mag2.madThresh);
    vidThreshYLim  = 3*std(vid.madThresh);
    threshYLabel = '|Velocity - median(Velocity)| (V/s)';
else
    error('"'+app.DesaccadeMethodDropDown.Value+'" is an invalid desaccading option.');
end

mag1VelThresh = mag1VelThreshSac;
mag2VelThresh = mag2VelThreshSac;
vidVelThresh  = vidVelThreshSac;

mag1VelThresh(mag1.saccades | mag1.saccades_manual) = nan;
mag2VelThresh(mag2.saccades | mag2.saccades_manual) = nan;
vidVelThresh(vid.saccades | vid.saccades_manual)    = nan;

mag1Pos = mag1.pos_data_aligned;
mag2Pos = mag2.pos_data_aligned;
vidPos  = vid.pos_data_upsampled_aligned;

mag1Pos(mag1.saccades | mag1.saccades_manual) = nan;
mag2Pos(mag2.saccades | mag2.saccades_manual) = nan;
vidPos(vid.saccades | vid.saccades_manual)    = nan;

mag1Vel = mag1.vel_data_aligned;
mag2Vel = mag2.vel_data_aligned;
vidVel  = vid.vel_data_upsampled_aligned;

mag1Vel(mag1.saccades | mag1.saccades_manual) = nan;
mag2Vel(mag2.saccades | mag2.saccades_manual) = nan;
vidVel(vid.saccades | vid.saccades_manual)    = nan;

mag1VelYLim = 6*mag1.vel_amp;
mag2VelYLim = 6*mag2.vel_amp;
vidVelYLim  = 6*vid.vel_amp; 

mag1Text = {sprintf('Fit Amp = %.3f deg/s', mag1.vel_amp), ...
                sprintf('Fit r^2 = %.3f', mag1.vel_fitr2), ...
                sprintf('Saccades = %.3f', mean(mag1.saccades_all)), ...
                sprintf('Scale Factor = %.3f deg/V', mag1.vel_scale)};
mag2Text = {sprintf('Fit Amp = %.3f deg/s', mag2.vel_amp), ...
                sprintf('Fit r^2 = %.3f', mag2.vel_fitr2), ...
                sprintf('Saccades = %.3f', mean(mag2.saccades_all)), ...
                sprintf('Scale Factor = %.3f deg/V', mag2.vel_scale)};
vidText = {sprintf('Fit Amp = %.3f deg/s', vid.vel_amp), ...
               sprintf('Fit r^2 = %.3f', vid.vel_fitr2), ...
               sprintf('Saccades = %.3f', mean(vid.saccades_all))};

sigColor = "#023dfa"; %"#1069bb"; 
sacColor = "#ff0511"; %[1 0.25 0.25]; 
altColor = "k";

if ~updateOnly
    %% Magnet 1 Threshold, Position, and Velocity
    cla(app.UIAxesMag1Threshold);
    cla(app.UIAxesMag1Position);
    cla(app.UIAxesMag1Velocity);
    
    timeLen = mag1.time_aligned(end);
    
    % Magnet 1 Threshold
    plot(app.UIAxesMag1Threshold, mag1.time_aligned, mag1VelThreshSac, 'Color',sacColor, 'DisplayName','Removed saccades');
    hold(app.UIAxesMag1Threshold, 'on');
    plot(app.UIAxesMag1Threshold, mag1.time_aligned, mag1VelThresh, 'Color',sigColor, 'DisplayName','Magnet 1');
    plot(app.UIAxesMag1Threshold, mag1.time_aligned, mag1Thresh, 'Color',altColor, 'LineStyle',':', 'LineWidth',2, 'DisplayName','Saccade Threshold');
    xlim(app.UIAxesMag1Threshold, [0, timeLen]);
    try
        ylim(app.UIAxesMag1Threshold, [0, mag1ThreshYLim]);
    catch
        % Do nothing
    end
    ylabel(app.UIAxesMag1Threshold, threshYLabel);
    legend(app.UIAxesMag1Threshold);
    hold(app.UIAxesMag1Threshold, 'off');
    
    % Magnet 1 Position
    plot(app.UIAxesMag1Position, mag1.time_aligned, mag1.pos_data_aligned, 'Color',sacColor, 'DisplayName','Removed saccades');
    hold(app.UIAxesMag1Position, 'on');
    plot(app.UIAxesMag1Position, mag1.time_aligned, mag1Pos, 'Color',sigColor, 'DisplayName','Magnet 1 Position');
    xlim(app.UIAxesMag1Position, [0, timeLen]);
    try
        ylim(app.UIAxesMag1Position, [min(mag1.pos_data_aligned), max(mag1.pos_data_aligned)]);
    catch
        % Do nothing
    end
    ylabel(app.UIAxesMag1Position, 'Position (V)');
    hold(app.UIAxesMag1Position, 'off');
    
    % Magnet 1 Velocity  
    plot(app.UIAxesMag1Velocity, mag1.time_aligned, mag1.vel_data_aligned, 'Color',sacColor, 'DisplayName','Removed saccades');
    hold(app.UIAxesMag1Velocity, 'on');
    plot(app.UIAxesMag1Velocity, mag1.time_aligned, mag1Vel, 'Color',sigColor, 'DisplayName','Magnet 1 Velocity');
    plot(app.UIAxesMag1Velocity, mag1.time_aligned, mag1.vel_fit, 'Color',altColor, 'LineWidth',2, 'DisplayName','Magnet 1 Velocity Fit');
    xlim(app.UIAxesMag1Velocity, [0, timeLen]);
    try
        ylim(app.UIAxesMag1Velocity, [-mag1VelYLim, mag1VelYLim]);
    catch
        % Do nothing
    end
    ylabel(app.UIAxesMag1Velocity, 'Velocity (V/s)');
    title(app.UIAxesMag1Velocity, mag1Text, ...
          'Units','characters', ...
          'Position',[6, 1, 0], ...
          'HorizontalAlignment','left', ...
          'FontSize',12, ...
          'BackgroundColor','w', ...
          'EdgeColor','k');
    hold(app.UIAxesMag1Velocity, 'off');
    
    % Link all x-axes
    linkaxes([app.UIAxesMag1Threshold, app.UIAxesMag1Position, app.UIAxesMag1Velocity], 'x');

    %% Magnet 2 Threshold, Position, and Velocity
    cla(app.UIAxesMag2Threshold);
    cla(app.UIAxesMag2Position);
    cla(app.UIAxesMag2Velocity);

    % Magnet 2 Threshold
    plot(app.UIAxesMag2Threshold, mag2.time_aligned, mag2VelThreshSac, 'Color',sacColor, 'DisplayName','Removed saccades');
    hold(app.UIAxesMag2Threshold, 'on');
    plot(app.UIAxesMag2Threshold, mag2.time_aligned, mag2VelThresh, 'Color',sigColor, 'DisplayName','Magnet 2');
    plot(app.UIAxesMag2Threshold, mag2.time_aligned, mag2Thresh, 'Color',altColor, 'LineStyle',':', 'LineWidth',2, 'DisplayName','Saccade Threshold');
    xlim(app.UIAxesMag2Threshold, [0, timeLen]);
    try
        ylim(app.UIAxesMag2Threshold, [0, mag2ThreshYLim]);
    catch
        % Do nothing
    end
    ylabel(app.UIAxesMag2Threshold, threshYLabel);
    legend(app.UIAxesMag2Threshold);
    hold(app.UIAxesMag2Threshold, 'off');
    
    % Magnet 2 Position
    plot(app.UIAxesMag2Position, mag2.time_aligned, mag2.pos_data_aligned, 'Color',sacColor, 'DisplayName','Removed saccades');
    hold(app.UIAxesMag2Position, 'on');
    plot(app.UIAxesMag2Position, mag2.time_aligned, mag2Pos, 'Color',sigColor, 'DisplayName','Magnet 2 Position');
    xlim(app.UIAxesMag2Position, [0, timeLen]);
    try
        ylim(app.UIAxesMag2Position, [min(mag2.pos_data_aligned), max(mag2.pos_data_aligned)]);
    catch
        % Do nothing
    end
    ylabel(app.UIAxesMag2Position, 'Position (V)');
    hold(app.UIAxesMag2Position, 'off');
    
    % Magnet 2 Velocity
    plot(app.UIAxesMag2Velocity, mag2.time_aligned, mag2.vel_data_aligned, 'Color',sacColor, 'DisplayName','Removed saccades');
    hold(app.UIAxesMag2Velocity, 'on');
    plot(app.UIAxesMag2Velocity, mag2.time_aligned, mag2Vel, 'Color',sigColor, 'DisplayName','Magnet 2 Velocity');
    plot(app.UIAxesMag2Velocity, mag2.time_aligned, mag2.vel_fit, 'Color',altColor, 'LineWidth',2, 'DisplayName','Magnet 2 Velocity Fit');
    xlim(app.UIAxesMag2Velocity, [0, timeLen]);
    ylim(app.UIAxesMag2Velocity, [-mag2VelYLim, mag2VelYLim]);
    ylabel(app.UIAxesMag2Velocity, 'Velocity (V/s)');
    title(app.UIAxesMag2Velocity, mag2Text, ...
          'Units','characters', ...
          'Position',[6, 1, 0], ...
          'HorizontalAlignment','left', ...
          'FontSize',12, ...
          'BackgroundColor','w', ...
          'EdgeColor','k');
    hold(app.UIAxesMag2Velocity, 'off');
    
    % Link all x-axes
    linkaxes([app.UIAxesMag2Threshold, app.UIAxesMag2Position, app.UIAxesMag2Velocity], 'x');

    %% Video Threshold, Position, and Velocity
    cla(app.UIAxesVidThreshold);
    cla(app.UIAxesVidPosition);
    cla(app.UIAxesVidVelocity);

    % Video Threshold
    plot(app.UIAxesVidThreshold, vid.time_upsampled_aligned, vidVelThreshSac, 'Color',sacColor, 'DisplayName','Removed saccades');
    hold(app.UIAxesVidThreshold, 'on');
    plot(app.UIAxesVidThreshold, vid.time_upsampled_aligned, vidVelThresh, 'Color',sigColor, 'DisplayName','Video Squared Velocity');
    plot(app.UIAxesVidThreshold, vid.time_upsampled_aligned, vidThresh, 'Color',altColor, 'LineStyle',':', 'LineWidth',2, 'DisplayName','Saccade Threshold');
    xlim(app.UIAxesVidThreshold, [0, timeLen]);
    ylim(app.UIAxesVidThreshold, [0, vidThreshYLim]);
    ylabel(app.UIAxesVidThreshold, threshYLabel);
    legend(app.UIAxesVidThreshold);
    hold(app.UIAxesVidThreshold, 'off');
    
    % Video Position
    plot(app.UIAxesVidPosition, vid.time_upsampled_aligned, vid.pos_data_upsampled_aligned, 'Color',sacColor, 'DisplayName','Removed saccades');
    hold(app.UIAxesVidPosition, 'on');
    plot(app.UIAxesVidPosition, vid.time_upsampled_aligned, vidPos, 'Color',sigColor, 'DisplayName','Video Position');
    xlim(app.UIAxesVidPosition, [0, timeLen]);
    ylim(app.UIAxesVidPosition, [min(vid.pos_data_upsampled_aligned), max(vid.pos_data_upsampled_aligned)]);
    ylabel(app.UIAxesVidPosition, 'Position (deg)');
    hold(app.UIAxesVidPosition, 'off');
    
    % Video Velocity
    plot(app.UIAxesVidVelocity, vid.time_upsampled_aligned, vid.vel_data_upsampled_aligned, 'Color',sacColor, 'DisplayName','Removed saccades');
    hold(app.UIAxesVidVelocity, 'on');
    plot(app.UIAxesVidVelocity, vid.time_upsampled_aligned, vidVel, 'Color',sigColor, 'DisplayName','Video Velocity');
    plot(app.UIAxesVidVelocity, vid.time_upsampled_aligned, vid.vel_fit, 'Color',altColor, 'LineWidth',2, 'DisplayName','Video Velocity Fit');
    xlim(app.UIAxesVidVelocity, [0, timeLen]);
    ylim(app.UIAxesVidVelocity, [-vidVelYLim, vidVelYLim]);
    ylabel(app.UIAxesVidVelocity, 'Velocity (deg/s)');
    title(app.UIAxesVidVelocity, vidText, ...
          'Units','characters', ...
          'Position',[6, 1, 0], ...
          'HorizontalAlignment','left', ...
          'FontSize',12, ...
          'BackgroundColor','w', ...
          'EdgeColor','k');
    hold(app.UIAxesVidVelocity, 'off');
    
    % Link all x-axes
    linkaxes([app.UIAxesVidThreshold, app.UIAxesVidPosition, app.UIAxesVidVelocity], 'x');
    
    %% Save data
    saveAnalysisInfo_APP;
else
    set(app.UIAxesMag1Threshold.Children(1), 'YData',mag1Thresh);
    set(app.UIAxesMag1Threshold.Children(2), 'YData',mag1VelThresh);
    set(app.UIAxesMag1Threshold.Children(3), 'YData',mag1VelThreshSac);
    set(app.UIAxesMag1Position.Children(1), 'YData',mag1Pos);
    set(app.UIAxesMag1Velocity.Children(2), 'YData',mag1Vel);
    try
        ylim(app.UIAxesMag1Threshold, [0, mag1ThreshYLim]);
        ylim(app.UIAxesMag1Velocity, [-mag1VelYLim, mag1VelYLim]);
    catch
        % Do nothing
    end
    app.UIAxesMag1Velocity.Title.String = mag1Text;

    set(app.UIAxesMag2Threshold.Children(1), 'YData',mag2Thresh);
    set(app.UIAxesMag2Threshold.Children(2), 'YData',mag2VelThresh);
    set(app.UIAxesMag2Threshold.Children(3), 'YData',mag2VelThreshSac);
    set(app.UIAxesMag2Position.Children(1), 'YData',mag2Pos);
    set(app.UIAxesMag2Velocity.Children(2), 'YData',mag2Vel);
    try
        ylim(app.UIAxesMag2Threshold, [0, mag2ThreshYLim]);
        ylim(app.UIAxesMag2Velocity, [-mag2VelYLim, mag2VelYLim]);
    catch
        % Do nothing
    end
    app.UIAxesMag2Velocity.Title.String = mag2Text;

    set(app.UIAxesVidThreshold.Children(1), 'YData',vidThresh);
    set(app.UIAxesVidThreshold.Children(2), 'YData',vidVelThresh);
    set(app.UIAxesVidThreshold.Children(3), 'YData',vidVelThreshSac);
    ylim(app.UIAxesVidThreshold, [0, vidThreshYLim]);
    set(app.UIAxesVidPosition.Children(1), 'YData',vidPos);
    set(app.UIAxesVidVelocity.Children(2), 'YData',vidVel);
    ylim(app.UIAxesVidVelocity, [-vidVelYLim, vidVelYLim]);
    app.UIAxesVidVelocity.Title.String = vidText;
end

end