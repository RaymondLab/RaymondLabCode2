function [vars] = saveplotSummaryDesaccading_APP(app, vars)
%SAVEPLOTSUMMARYDESACCADING_APP

%% Load data
loadAnalysisInfo_APP;

mag1Vel = mag1.vel_data_aligned;
mag2Vel = mag2.vel_data_aligned;
vidVel  = vid.vel_data_upsampled_aligned;

mag1Vel(mag1.saccades_all) = nan;
mag2Vel(mag2.saccades_all) = nan;
vidVel(vid.saccades_all)   = nan;

set(app.UIAxesMag1Velocity.Children(2), 'YData',mag1Vel);
set(app.UIAxesMag2Velocity.Children(2), 'YData',mag2Vel);
set(app.UIAxesVidVelocity.Children(2), 'YData',vidVel);

%% Figure: Desaccaded Position and Velocity Traces
fig = figure('units','normalized', 'outerposition',[0 0 1 1]);
h = tiledlayout(fig, 3, 1, 'TileSpacing','tight', 'Padding','compact');
ax1 = nexttile;
copyobj(app.UIAxesMag1Velocity.Children, ax1);
legend(ax1, {'Removed Saccades', 'Velocity Trace', 'Sinusoidal Fit'});
ax1.XLim = app.UIAxesMag1Velocity.XLim;
ax1.YLim = app.UIAxesMag1Velocity.YLim;
ax1.YLabel.String = app.UIAxesMag1Velocity.YLabel.String;
title(ax1, 'Magnet Channel 1: Unscaled Velocity', 'FontSize',16);
mag1Text = {sprintf('Fit Amp = %.3f deg/s', mag1.vel_amp), ...
            sprintf('Fit r^2 = %.3f', mag1.vel_fitr2), ...
            sprintf('Saccades = %.3f', mean(mag1.saccades_all)), ...
            sprintf('Scale Factor = %.3f deg/V', mag1.vel_scale)};
annotation(fig, 'textbox',[0.045 0.656 0.3 0.3], 'String',mag1Text, 'FitBoxToText','on', 'BackgroundColor','w');
grid(ax1, 'on');
box(ax1, 'on');

ax2 = nexttile;
copyobj(app.UIAxesVidVelocity.Children, ax2);
ax2.XLim = app.UIAxesVidVelocity.XLim;
ax2.YLim = app.UIAxesVidVelocity.YLim;
ax2.YLabel.String = app.UIAxesVidVelocity.YLabel.String;
title(ax2, 'Video Channel: Desaccaded Velocity', 'FontSize',16);
vidText = {sprintf('Fit Amp = %.3f deg/s', vid.vel_amp), ...
           sprintf('Fit r^2 = %.3f', vid.vel_fitr2), ...
           sprintf('Saccades = %.3f', mean(vid.saccades_all))};
annotation(fig, 'textbox',[0.045 0.34 0.3 0.3], 'String',vidText, 'FitBoxToText','on', 'BackgroundColor','w');
grid(ax2, 'on');
box(ax2, 'on');

ax3 = nexttile;
copyobj(app.UIAxesMag2Velocity.Children, ax3);
ax3.XLim = app.UIAxesMag2Velocity.XLim;
ax3.YLim = app.UIAxesMag2Velocity.YLim;
ax3.XLabel.String = 'Time (s)';
ax3.YLabel.String = app.UIAxesMag2Velocity.YLabel.String;
title(ax3, 'Magnet Channel 2: Unscaled Velocity', 'FontSize',16);
mag2Text = {sprintf('Fit Amp = %.3f deg/s', mag2.vel_amp), ...
            sprintf('Fit r^2 = %.3f', mag2.vel_fitr2), ...
            sprintf('Saccades = %.3f', mean(mag2.saccades_all)), ...
            sprintf('Scale Factor = %.3f deg/V', mag2.vel_scale)};
annotation(fig, 'textbox',[0.045 0.025 0.3 0.3], 'String',mag2Text, 'FitBoxToText','on', 'BackgroundColor','w');
grid(ax3, 'on');
box(ax3, 'on');

% Link the x-axes
linkaxes([ax1, ax2, ax3], 'x');

%% Save Figures
set(fig, 'Units','Inches');
pos = get(fig, 'Position');
set(fig,'PaperPositionMode','Auto', 'PaperUnits','Inches', 'PaperSize',[pos(3),pos(4)]);
print(fig, fullfile(cd,'Summary_Desaccading.pdf'), '-vector', '-dpdf');
savefig('Summary_Desaccading.fig');

%% Save data
saveAnalysisInfo_APP;

pause(0.01);
close(fig);

end