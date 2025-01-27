function vars = saveplotLinearityAnalysis_APP(app, vars)
%SAVEPLOTLINEARITYANALYSIS_APP

%% Load data
loadAnalysisInfo_APP;

%% Linearity Calculation
c = linspace(1, 10, length(mag1.pos_data_aligned_scaledInVel));
vidTime = vid.time_upsampled_aligned;
mag1Time = mag1.time_aligned;
mag2Time = mag2.time_aligned;

% Get saccade start and end times, exit early upon error
[mag1GoodStarts,mag1GoodStops,errInfo1] = getSaccadeStartEndTimes(mag1.saccades_all);
[mag2GoodStarts,mag2GoodStops,errInfo2] = getSaccadeStartEndTimes(mag2.saccades_all);
if errInfo1 || errInfo2
    return
end

% Get aligned desaccaded position data
mag1_pos_aligned = mag1.pos_data_aligned_scaledInVel;
mag1_pos_aligned(mag1.saccades_all) = nan;
vid1_pos_aligned = vid.pos_data_upsampled_aligned;
vid1_pos_aligned(mag1.saccades_all) = nan;

mag2_pos_aligned = mag2.pos_data_aligned_scaledInVel;
mag2_pos_aligned(mag2.saccades_all) = nan;
vid2_pos_aligned = vid.pos_data_upsampled_aligned;
vid2_pos_aligned(mag2.saccades_all) = nan;

% Get aligned desaccaded velocity data
mag1_vel_aligned = mag1.vel_data_aligned_scaledInVel;
mag1_vel_aligned(mag1.saccades_all) = nan;
vid1_vel_aligned = vid.vel_data_upsampled_aligned;
vid1_vel_aligned(mag1.saccades_all) = nan;

mag2_vel_aligned = mag2.vel_data_aligned_scaledInVel;
mag2_vel_aligned(mag2.saccades_all) = nan;
vid2_vel_aligned = vid.vel_data_upsampled_aligned;
vid2_vel_aligned(mag2.saccades_all) = nan;

% Measure linearity of whole segments
mag1.pos_lin_all = measureLinearity(mag1.pos_data_aligned_scaledInVel(~mag1.saccades_all), vid.pos_data_upsampled_aligned(~mag1.saccades_all), ~mag1.saccades_all);
mag1.vel_lin_all = measureLinearity(mag1.vel_data_aligned_scaledInVel(~mag1.saccades_all), vid.vel_data_upsampled_aligned(~mag1.saccades_all), ~mag1.saccades_all);
mag2.pos_lin_all = measureLinearity(mag2.pos_data_aligned_scaledInVel(~mag2.saccades_all), vid.pos_data_upsampled_aligned(~mag2.saccades_all), ~mag2.saccades_all);
mag2.vel_lin_all = measureLinearity(mag2.vel_data_aligned_scaledInVel(~mag2.saccades_all), vid.vel_data_upsampled_aligned(~mag2.saccades_all), ~mag2.saccades_all);

% Measure linearity of non-saccades, split up into chunks
for i = 1:length(mag1GoodStarts)
    chunk = mag1GoodStarts(i):mag1GoodStops(i);
    mag1.pos_lin_chunks(i) = measureLinearity(mag1.pos_data_aligned_scaledInVel(chunk), vid.pos_data_upsampled_aligned(chunk));
    mag1.vel_lin_chunks(i) = measureLinearity(mag1.vel_data_aligned_scaledInVel(chunk), vid.vel_data_upsampled_aligned(chunk));
end

for i = 1:length(mag2GoodStarts)
    chunk = mag2GoodStarts(i):mag2GoodStops(i);
    mag2.pos_lin_chunks(i) = measureLinearity(mag2.pos_data_aligned_scaledInVel(chunk), vid.pos_data_upsampled_aligned(chunk));
    mag2.vel_lin_chunks(i) = measureLinearity(mag2.vel_data_aligned_scaledInVel(chunk), vid.vel_data_upsampled_aligned(chunk));
end


%% Figure 1: Linearity of Position Traces of both Magnet Channels
fig1 = figure('units','normalized', 'outerposition',[0 0 1 1]); clf
h1 = tiledlayout(fig1, 2, 3, 'TileSpacing','tight', 'Padding','compact');

% Magnet Channel 1, Fit Coefficients
ax1 = nexttile;
for i = 1:length(mag1GoodStarts)
    chunk = mag1GoodStarts(i):mag1GoodStops(i);

    magPoints = mag1_pos_aligned(chunk);
    magPoints = magPoints - nanmean(magPoints);
    mag1_pos_aligned(chunk) = magPoints;

    vidChunk = vid1_pos_aligned(chunk);
    vidChunk = vidChunk - nanmean(vidChunk);
    vid1_pos_aligned(chunk) = vidChunk;
end

scatter(ax1, mag1_pos_aligned, vid1_pos_aligned, 50, c, '.', 'HandleVisibility','off');
hold(ax1, 'on');
xylimit = max(max(abs(mag1_pos_aligned)), max(abs(vid1_pos_aligned)));

coefficients = polyfit(mag1_pos_aligned(~mag1.saccades_all), vid1_pos_aligned(~mag1.saccades_all), 1);
slopes = [];
for i = 1:length(mag1GoodStarts)
    chunk = mag1GoodStarts(i)+1:mag1GoodStops(i);

    magPoints = mag1_pos_aligned(chunk);
    magPoints = magPoints - nanmean(magPoints);

    vidChunk = vid1_pos_aligned(chunk);
    vidChunk = vidChunk - nanmean(vidChunk);
    coeffs = polyfit(magPoints, vidChunk, 1);
    slopes(end+1) = coeffs(1) - coefficients(1);
    xFit = linspace(-xylimit, xylimit, 1000);
    yFit = polyval(coeffs , xFit);
    plot(ax1, xFit, yFit, '-k', 'LineWidth', .1, 'HandleVisibility','off');
end
mstd = std(slopes);
mvar = var(slopes);

xFit = linspace(-xylimit, xylimit, 1000);
yFit = polyval(coefficients , xFit);
plot(ax1, xFit, yFit, '-r', 'LineWidth',3, 'DisplayName',sprintf(' m = %.5f (std=%.3f, var=%.3f)', coefficients(1), mstd, mvar));
xline(ax1, 0, ':r', 'HandleVisibility','off');
yline(ax1, 0, ':r', 'HandleVisibility','off');
ylim(ax1, [-xylimit, xylimit]);
xlim(ax1, [-xylimit, xylimit]);
ax1.XLabel.String = 'Magnet Channel 1 Position';
ax1.YLabel.String = 'Video Channel Position';
title(ax1, sprintf('Magnet Channel 1: (%d) Non-Saccade Segments',length(slopes)), 'FontSize',14);
grid(ax1, 'on');
box(ax1, 'on');
hold(ax1, 'off');
legend(ax1, 'Location','northwest', 'FontSize',12);

ax2 = nexttile;
allData = logical(ones(length(~mag1.saccades_all), 1));
hold(ax2, 'on');
Rsq = linearityScatterPlot(ax2, mag1.pos_data_aligned_scaledInVel, vid.pos_data_upsampled_aligned, allData, c);
[magmin, magmax] = bounds(mag1.pos_data_aligned_scaledInVel(allData));
[vidmin, vidmax] = bounds(vid.pos_data_upsampled_aligned(allData));
title(ax2, sprintf('Magnet Channel 1 (with saccades): r^2 = %.5f', Rsq), 'FontSize',14);
xlabel(ax2, sprintf('Magnet Channel 1 Position | minmax = (%.1f, %.1f)', magmin, magmax));
ylabel(ax2, sprintf('Video Position | minmax = (%.1f, %.1f)', vidmin, vidmax));
grid(ax2, 'on');
box(ax2, 'on');
hold(ax2, 'off');

ax3 = nexttile;
scatter(ax3, mag1.pos_data_aligned_scaledInVel(~mag1.saccades_all), vid.pos_data_upsampled_aligned(~mag1.saccades_all), 4, c(~mag1.saccades_all), '.');
[magmin, magmax] = bounds(mag1.pos_data_aligned_scaledInVel(~mag1.saccades_all));
[vidmin, vidmax] = bounds(vid.pos_data_upsampled_aligned(~mag1.saccades_all));
llimits = min(magmin, vidmin);
ulimits = max(magmax, vidmax);
hold(ax3, 'on');
plot(ax3, mag1.pos_lin_all.range, mag1.pos_lin_all.yfit, 'k', 'lineWidth', 2);
colormap(ax3, hsv);
title(ax3, sprintf('Magnet Channel 1 (saccades removed): r^2 = %.5f', mag1.pos_lin_all.Rsq), 'FontSize',14)
xlabel(ax3, sprintf('Magnet Channel 1 Position | minmax = (%.1f, %.1f)', magmin, magmax));
ylabel(ax3, sprintf('Video Position | minmax = (%.1f, %.1f)', vidmin, vidmax));
grid(ax3, 'on');
box(ax3, 'on');
ylim(ax3, [llimits, ulimits]);
xlim(ax3, [llimits, ulimits]);
hold(ax3, 'off');

% Magnet Channel 2, Fit Coefficients
ax4 = nexttile;
for i = 1:length(mag2GoodStarts)
    chunk = mag2GoodStarts(i):mag2GoodStops(i);

    magPoints = mag2_pos_aligned(chunk);
    magPoints = magPoints - nanmean(magPoints);
    mag2_pos_aligned(chunk) = magPoints;

    vidChunk = vid2_pos_aligned(chunk);
    vidChunk = vidChunk - nanmean(vidChunk);
    vid2_pos_aligned(chunk) = vidChunk;
end

scatter(ax4, mag2_pos_aligned, vid2_pos_aligned, 50, c, '.', 'HandleVisibility','off');
hold(ax4, 'on');
xylimit = max(max(abs(mag2_pos_aligned)), max(abs(vid2_pos_aligned)));

coefficients = polyfit(mag2_pos_aligned(~mag2.saccades_all), vid2_pos_aligned(~mag2.saccades_all), 1);
slopes = [];
for i = 1:length(mag2GoodStarts)
    chunk = mag2GoodStarts(i)+1:mag2GoodStops(i);

    magPoints = mag2_pos_aligned(chunk);
    magPoints = magPoints - nanmean(magPoints);

    vidChunk = vid2_pos_aligned(chunk);
    vidChunk = vidChunk - nanmean(vidChunk);
    coeffs = polyfit(magPoints, vidChunk, 1);
    slopes(end+1) = coeffs(1) - coefficients(1);
    xFit = linspace(-xylimit, xylimit, 1000);
    yFit = polyval(coeffs , xFit);
    plot(ax4, xFit, yFit, '-k', 'LineWidth', .1, 'HandleVisibility','off');
end
mstd = nanstd(slopes);
mvar = nanvar(slopes);

xFit = linspace(-xylimit, xylimit, 1000);
yFit = polyval(coefficients , xFit);
plot(ax4, xFit, yFit, '-r', 'LineWidth', 3, 'DisplayName',sprintf(' m = %.5f (std=%.3f, var=%.3f)', coefficients(1), mstd, mvar));
xline(ax4, 0, ':r', 'HandleVisibility','off');
yline(ax4, 0, ':r', 'HandleVisibility','off');
ylim(ax4, [-xylimit, xylimit]);
xlim(ax4, [-xylimit, xylimit]);
ax4.XLabel.String = 'Magnet Channel 2 Position';
ax4.YLabel.String = 'Video Channel Position';
title(ax4, sprintf('Magnet Channel 2: (%d) Non-Saccade Segments',length(slopes)), 'FontSize',14);
grid(ax4, 'on');
box(ax4, 'on');
hold(ax4, 'off');
legend(ax4, 'Location','northwest', 'FontSize',12);

ax5 = nexttile;
allData = logical(ones(length(~mag2.saccades_all), 1));
hold(ax5, 'on');
Rsq = linearityScatterPlot(ax5, mag2.pos_data_aligned_scaledInVel, vid.pos_data_upsampled_aligned, allData, c);
[magmin, magmax] = bounds(mag2.pos_data_aligned_scaledInVel(allData));
[vidmin, vidmax] = bounds(vid.pos_data_upsampled_aligned(allData));
title(ax5, sprintf('Magnet Channel 2 (with saccades): r^2 = %.5f', Rsq), 'FontSize',14);
xlabel(ax5, sprintf('Magnet Channel 2 Position | minmax = (%.1f, %.1f)', magmin, magmax));
ylabel(ax5, sprintf('Video Position | minmax = (%.1f, %.1f)', vidmin, vidmax));
grid(ax5, 'on');
box(ax5, 'on');
hold(ax5, 'off');

ax6 = nexttile;
scatter(ax6, mag2.pos_data_aligned_scaledInVel(~mag2.saccades_all), vid.pos_data_upsampled_aligned(~mag2.saccades_all), 4, c(~mag2.saccades_all), '.');
[magmin, magmax] = bounds(mag2.pos_data_aligned_scaledInVel(~mag2.saccades_all));
[vidmin, vidmax] = bounds(vid.pos_data_upsampled_aligned(~mag2.saccades_all));
llimits = min(magmin, vidmin);
ulimits = max(magmax, vidmax);
hold(ax6, 'on');
plot(ax6, mag2.pos_lin_all.range, mag2.pos_lin_all.yfit, 'k', 'lineWidth', 2);
colormap(ax6, hsv);
title(ax6, sprintf('Magnet Channel 2 (saccades removed): r^2 = %.5f', mag2.pos_lin_all.Rsq), 'FontSize',14)
xlabel(ax6, sprintf('Magnet Channel 2 Position | minmax = (%.1f, %.1f)', magmin, magmax));
ylabel(ax6, sprintf('Video Position | minmax = (%.1f, %.1f)', vidmin, vidmax));
grid(ax6, 'on');
box(ax6, 'on');
ylim(ax6, [llimits, ulimits]);
xlim(ax6, [llimits, ulimits]);
hold(ax6, 'off');

title(h1, 'Linearity (Rainbow) Plots', 'FontSize',18);

%% Save Figure
set(fig1, 'Units','Inches');
pos = get(fig1, 'Position');
set(fig1,'PaperPositionMode','Auto', 'PaperUnits','Inches', 'PaperSize',[pos(3),pos(4)]);
print(fig1, fullfile(cd,'Summary_LinearityAnalysis_Main_Position.pdf'), '-vector', '-dpdf');
savefig(fig1, 'Summary_LinearityAnalysis_Main_Position.fig');


% %% Figure 2: Linearity of Velocity Traces of both Magnet Channels
% fig2 = figure('units','normalized', 'outerposition',[0 0 1 1]); clf
% h2 = tiledlayout(fig2, 2, 3, 'TileSpacing','tight', 'Padding','compact');
% 
% % Magnet Channel 1, Fit Coefficients
% ax7 = nexttile;
% for i = 1:length(mag1GoodStarts)
%     chunk = mag1GoodStarts(i):mag1GoodStops(i);
% 
%     magPoints = mag1_vel_aligned(chunk);
%     magPoints = magPoints - nanmean(magPoints);
%     mag1_vel_aligned(chunk) = magPoints;
% 
%     vidChunk = vid1_vel_aligned(chunk);
%     vidChunk = vidChunk - nanmean(vidChunk);
%     vid1_vel_aligned(chunk) = vidChunk;
% end
% 
% scatter(ax7, mag1_vel_aligned, vid1_vel_aligned, 50, c, '.', 'HandleVisibility','off');
% hold(ax7, 'on');
% xylimit = max(max(abs(mag1_vel_aligned)), max(abs(vid1_vel_aligned)));
% 
% coefficients = polyfit(mag1_vel_aligned(~mag1.saccades_all), vid1_vel_aligned(~mag1.saccades_all), 1);
% slopes = [];
% for i = 1:length(mag1GoodStarts)
%     chunk = mag1GoodStarts(i)+1:mag1GoodStops(i);
% 
%     magPoints = mag1_vel_aligned(chunk);
%     magPoints = magPoints - nanmean(magPoints);
% 
%     vidChunk = vid1_vel_aligned(chunk);
%     vidChunk = vidChunk - nanmean(vidChunk);
%     coeffs = polyfit(magPoints, vidChunk, 1);
%     slopes(end+1) = coeffs(1) - coefficients(1);
%     xFit = linspace(-xylimit, xylimit, 1000);
%     yFit = polyval(coeffs , xFit);
%     plot(ax7, xFit, yFit, '-k', 'LineWidth', .1, 'HandleVisibility','off');
% end
% mstd = std(slopes);
% mvar = var(slopes);
% 
% xFit = linspace(-xylimit, xylimit, 1000);
% yFit = polyval(coefficients , xFit);
% plot(ax7, xFit, yFit, '-r', 'LineWidth',3, 'DisplayName',sprintf(' m = %.5f (std=%.3f, var=%.3f)', coefficients(1), mstd, mvar));
% xline(ax7, 0, ':r', 'HandleVisibility','off');
% yline(ax7, 0, ':r', 'HandleVisibility','off');
% ylim(ax7, [-xylimit, xylimit]);
% xlim(ax7, [-xylimit, xylimit]);
% ax7.XLabel.String = 'Magnet Channel 1 Velocity';
% ax7.YLabel.String = 'Video Channel Velocity';
% title(ax7, sprintf('Magnet Channel 1: (%d) Non-Saccade Segments',length(slopes)), 'FontSize',14);
% grid(ax7, 'on');
% box(ax7, 'on');
% hold(ax7, 'off');
% legend(ax7, 'Location','northwest', 'FontSize',12);
% 
% ax8 = nexttile;
% allData = logical(ones(length(~mag1.saccades_all), 1));
% hold(ax8, 'on');
% Rsq = linearityScatterPlot(ax8, mag1.vel_data_aligned_scaledInVel, vid.vel_data_upsampled_aligned, allData, c);
% [magmin, magmax] = bounds(mag1.vel_data_aligned_scaledInVel(allData));
% [vidmin, vidmax] = bounds(vid.vel_data_upsampled_aligned(allData));
% title(ax8, ['Magnet Channel 1 (with saccades): r^2: ', num2str(Rsq)], 'FontSize',14);
% xlabel(ax8, sprintf('Magnet Channel 1 Velocity | minmax = (%.1f, %.1f)', magmin, magmax));
% ylabel(ax8, sprintf('Video Velocity | minmax = (%.1f, %.1f)', vidmin, vidmax));
% grid(ax8, 'on');
% box(ax8, 'on');
% hold(ax8, 'off');
% 
% ax9 = nexttile;
% scatter(ax9, mag1.vel_data_aligned_scaledInVel(~mag1.saccades_all), vid.vel_data_upsampled_aligned(~mag1.saccades_all), 4, c(~mag1.saccades_all), '.');
% [magmin, magmax] = bounds(mag1.vel_data_aligned_scaledInVel(~mag1.saccades_all));
% [vidmin, vidmax] = bounds(vid.vel_data_upsampled_aligned(~mag1.saccades_all));
% hold(ax9, 'on');
% plot(ax9, mag1.vel_lin_all.range, mag1.vel_lin_all.yfit, 'k', 'lineWidth', 2);
% colormap(ax9, hsv);
% title(ax9, ['Magnet Channel 1 (saccades removed): r^2: ', num2str(mag1.vel_lin_all.Rsq)], 'FontSize',14)
% xlabel(ax9, sprintf('Magnet Channel 1 Velocity | minmax = (%.1f, %.1f)', magmin, magmax));
% ylabel(ax9, sprintf('Video Velocity | minmax = (%.1f, %.1f)', vidmin, vidmax));
% grid(ax9, 'on');
% box(ax9, 'on');
% hold(ax9, 'off');
% llimits = min(magmin, vidmin);
% ulimits = max(magmax, vidmax);
% ylim(ax9, [llimits, ulimits]);
% xlim(ax9, [llimits, ulimits]);
% 
% % Magnet Channel 2, Fit Coefficients
% ax10 = nexttile;
% for i = 1:length(mag2GoodStarts)
%     chunk = mag2GoodStarts(i):mag2GoodStops(i);
% 
%     magPoints = mag2_vel_aligned(chunk);
%     magPoints = magPoints - nanmean(magPoints);
%     mag2_vel_aligned(chunk) = magPoints;
% 
%     vidChunk = vid2_vel_aligned(chunk);
%     vidChunk = vidChunk - nanmean(vidChunk);
%     vid2_vel_aligned(chunk) = vidChunk;
% end
% 
% scatter(ax10, mag2_vel_aligned, vid2_vel_aligned, 50, c, '.', 'HandleVisibility','off');
% hold(ax10, 'on');
% xylimit = max(max(abs(mag2_vel_aligned)), max(abs(vid2_vel_aligned)));
% 
% coefficients = polyfit(mag2_vel_aligned(~mag2.saccades_all), vid2_vel_aligned(~mag2.saccades_all), 1);
% slopes = [];
% for i = 1:length(mag2GoodStarts)
%     chunk = mag2GoodStarts(i)+1:mag2GoodStops(i);
% 
%     magPoints = mag2_vel_aligned(chunk);
%     magPoints = magPoints - nanmean(magPoints);
% 
%     vidChunk = vid2_vel_aligned(chunk);
%     vidChunk = vidChunk - nanmean(vidChunk);
%     coeffs = polyfit(magPoints, vidChunk, 1);
%     slopes(end+1) = coeffs(1) - coefficients(1);
%     xFit = linspace(-xylimit, xylimit, 1000);
%     yFit = polyval(coeffs , xFit);
%     plot(ax10, xFit, yFit, '-k', 'LineWidth', .1, 'HandleVisibility','off');
% end
% mstd = nanstd(slopes);
% mvar = nanvar(slopes);
% 
% xFit = linspace(-xylimit, xylimit, 1000);
% yFit = polyval(coefficients , xFit);
% plot(ax10, xFit, yFit, '-r', 'LineWidth', 3, 'DisplayName',sprintf(' m = %.5f (std=%.3f, var=%.3f)', coefficients(1), mstd, mvar));
% xline(ax10, 0, ':r', 'HandleVisibility','off');
% yline(ax10, 0, ':r', 'HandleVisibility','off');
% ylim(ax10, [-xylimit, xylimit]);
% xlim(ax10, [-xylimit, xylimit]);
% ax10.XLabel.String = 'Magnet Channel 2 Velocity';
% ax10.YLabel.String = 'Video Channel Velocity';
% title(ax10, sprintf('Magnet Channel 2: (%d) Non-Saccade Segments',length(slopes)), 'FontSize',14);
% grid(ax10, 'on');
% box(ax10, 'on');
% hold(ax10, 'off');
% legend(ax10, 'Location','northwest', 'FontSize',12);
% 
% ax11 = nexttile;
% allData = logical(ones(length(~mag2.saccades_all), 1));
% hold(ax11, 'on');
% Rsq = linearityScatterPlot(ax11, mag2.vel_data_aligned_scaledInVel, vid.vel_data_upsampled_aligned, allData, c);
% [magmin, magmax] = bounds(mag2.vel_data_aligned_scaledInVel(allData));
% [vidmin, vidmax] = bounds(vid.vel_data_upsampled_aligned(allData));
% title(ax11, ['Magnet Channel 2 (with saccades): r^2: ', num2str(Rsq)], 'FontSize',14);
% xlabel(ax11, sprintf('Magnet Channel 2 Velocity | minmax = (%.1f, %.1f)', magmin, magmax));
% ylabel(ax11, sprintf('Video Velocity | minmax = (%.1f, %.1f)', vidmin, vidmax));
% grid(ax11, 'on');
% box(ax11, 'on');
% hold(ax11, 'off');
% 
% ax12 = nexttile;
% scatter(ax12, mag2.vel_data_aligned_scaledInVel(~mag2.saccades_all), vid.vel_data_upsampled_aligned(~mag2.saccades_all), 4, c(~mag2.saccades_all), '.');
% [magmin, magmax] = bounds(mag2.vel_data_aligned_scaledInVel(~mag2.saccades_all));
% [vidmin, vidmax] = bounds(vid.vel_data_upsampled_aligned(~mag2.saccades_all));
% hold(ax12, 'on');
% plot(ax12, mag2.vel_lin_all.range, mag2.vel_lin_all.yfit, 'k', 'lineWidth', 2);
% colormap(ax12, hsv);
% title(ax12, ['Magnet Channel 2 (saccades removed): r^2: ', num2str(mag2.vel_lin_all.Rsq)], 'FontSize',14)
% xlabel(ax12, sprintf('Magnet Channel 2 Velocity | minmax = (%.1f, %.1f)', magmin, magmax));
% ylabel(ax12, sprintf('Video Velocity | minmax = (%.1f, %.1f)', vidmin, vidmax));
% grid(ax12, 'on');
% box(ax12, 'on');
% hold(ax12, 'off');
% llimits = min(magmin, vidmin);
% ulimits = max(magmax, vidmax);
% ylim(ax12, [llimits, ulimits]);
% xlim(ax12, [llimits, ulimits]);
% 
% %% Save Figure
% set(fig2, 'Units','Inches');
% pos = get(fig2, 'Position');
% set(fig2,'PaperPositionMode','Auto', 'PaperUnits','Inches', 'PaperSize',[pos(3),pos(4)]);
% print(fig2, fullfile(cd,'Summary_LinearityAnalysis_Velocity.pdf'), '-vector', '-dpdf');
% savefig(fig2, 'Summary_LinearityAnalysis_Velocity.fig');


%% Figure 3: Bland Altman Plot of Position Traces
fig3 = figure('units','normalized', 'outerposition',[0 0 0.9 0.8]); clf
h3 = tiledlayout(fig3, 1, 2, 'TileSpacing','compact', 'Padding','compact');

xavgs1 = 0.5 * (mag1.pos_data_aligned_scaledInVel + vid.pos_data_upsampled_aligned);
xavgs1(mag1.saccades_all) = nan;
xavgs2 = 0.5 * (mag2.pos_data_aligned_scaledInVel + vid.pos_data_upsampled_aligned);
xavgs2(mag2.saccades_all) = nan;
ydiffs1 = mag1.pos_data_aligned_scaledInVel - vid.pos_data_upsampled_aligned;
ydiffs1(mag1.saccades_all) = nan;
ydiffs2 = mag1.pos_data_aligned_scaledInVel - vid.pos_data_upsampled_aligned;
ydiffs2(mag2.saccades_all) = nan;

% Prepare in advance plot boundaries. They will be needed afterwards
delta_x = abs(0.3 * max(xavgs1));
x_min = floor(min(xavgs1) - delta_x);
x_max = floor(max(xavgs1) + delta_x);
delta_y = abs(0.3 * max(ydiffs1));
y_min = floor(min(ydiffs1) - delta_y);
y_max = floor(max(ydiffs1) + delta_y);
xlims1 = [x_min, x_max];
ylims1 = [y_min, y_max];

delta_x = abs(0.3 * max(xavgs2));
x_min = floor(min(xavgs2) - delta_x);
x_max = floor(max(xavgs2) + delta_x);
delta_y = abs(0.3 * max(ydiffs2));
y_min = floor(min(ydiffs2) - delta_y);
y_max = floor(max(ydiffs2) + delta_y);
xlims2 = [x_min, x_max];
ylims2 = [y_min, y_max];

% Compute Bias and UB/LB
b1 = nanmean(ydiffs1);
UB1 = b1 + 1.96 * nanstd(ydiffs1);
LB1 = b1 - 1.96 * nanstd(ydiffs1);

b2 = nanmean(ydiffs2);
UB2 = b2 + 1.96 * nanstd(ydiffs2);
LB2 = b2 - 1.96 * nanstd(ydiffs2);

P1a = polyfit(xavgs1(~mag1.saccades_all), ydiffs1(~mag1.saccades_all), 1);
P1b = polyfit(xavgs1(~mag1.saccades_all), ydiffs1(~mag1.saccades_all), 3);
xfit1 = linspace(min(xavgs1), max(xavgs1), 1000);
yfit1a = polyval(P1a, xfit1);
yfit1b = polyval(P1b, xfit1);
txt1a = sprintf(' y = %.4f * x + %.4f', P1a(1), P1a(2));
txt1b = sprintf(' y = %.4f * x^3 + %.4f * x^2 + %.4f * x + %.4f', P1b(1), P1b(2), P1b(3), P1b(4));

P2a = polyfit(xavgs2(~mag2.saccades_all), ydiffs2(~mag2.saccades_all), 1);
P2b = polyfit(xavgs2(~mag2.saccades_all), ydiffs2(~mag2.saccades_all), 3);
xfit2 = linspace(min(xavgs2), max(xavgs2), 1000);
yfit2a = polyval(P2a, xfit2);
yfit2b = polyval(P2b, xfit2);
txt2a = sprintf(' y = %.4f * x + %.4f', P2a(1), P2a(2));
txt2b = sprintf(' y = %.4f * x^3 + %.4f * x^2 + %.4f * x + %.4f', P2b(1), P2b(2), P2b(3), P2b(4));


% Magnet Channel 1
ax13 = nexttile;
hold(ax13, 'on');
scatter(ax13, xavgs1(~mag1.saccades_all), ydiffs1(~mag1.saccades_all), 4, c(~mag1.saccades_all), '.', 'HandleVisibility','off');
yline(ax13, b1, '--k', sprintf('bias = %.4f',b1), 'LineWidth',1.5, 'FontSize',14, 'HandleVisibility','off');
yline(ax13, UB1, '--k', sprintf('UB = %.4f',UB1), 'LineWidth',1.5, 'FontSize',14, 'HandleVisibility','off');
yline(ax13, LB1, '--k', sprintf('LB = %.4f',LB1), 'LineWidth',1.5, 'FontSize',14, 'HandleVisibility','off');
plot(ax13, xfit1, yfit1a, '-k', 'LineWidth',3, 'DisplayName',txt1a);
plot(ax13, xfit1, yfit1b, '-r', 'LineWidth',2, 'DisplayName',txt1b);
xlabel(ax13, 'Mean of MAG1 and VIDEO positions ');
ylabel(ax13, 'MAG1 - VIDEO (position)');
xlim(ax13, xlims1);
ylim(ax13, ylims1);
hold(ax13, 'off');
grid(ax13, 'on');
box(ax13, 'on');
title(ax13, 'Magnet Channel 1 Position (saccades removed)', 'FontSize',14);
legend(ax13, 'Location','northwest', 'FontSize',12);

ax14 = nexttile;
hold(ax14, 'on');
scatter(ax14, xavgs2(~mag2.saccades_all), ydiffs2(~mag2.saccades_all), 4, c(~mag2.saccades_all), '.', 'HandleVisibility','off');
yline(ax14, b2, '--k', sprintf('bias = %.4f',b2), 'LineWidth',1.5, 'FontSize',14, 'HandleVisibility','off');
yline(ax14, UB2, '--k', sprintf('UB = %.4f',UB2), 'LineWidth',1.5, 'FontSize',14, 'HandleVisibility','off');
yline(ax14, LB2, '--k', sprintf('LB = %.4f',LB2), 'LineWidth',1.5, 'FontSize',14, 'HandleVisibility','off');
plot(ax14, xfit2, yfit2a, '-k', 'LineWidth',3, 'DisplayName',txt2a);
plot(ax14, xfit2, yfit2b, '-r', 'LineWidth',2, 'DisplayName',txt2b);
xlabel(ax14, 'Mean of MAG2 and VIDEO positions ');
ylabel(ax14, 'MAG2 - VIDEO (position)');
xlim(ax14, xlims2);
ylim(ax14, ylims2);
hold(ax14, 'off');
grid(ax14, 'on');
box(ax14, 'on');
title(ax14, 'Magnet Channel 2 Position (saccades removed)', 'FontSize',14);
legend(ax14, 'Location','northwest', 'FontSize',12);

title(h3, 'Bland-Altman Plots (Relative to Video Position)', 'FontSize',18);

% Save Figure
set(fig3, 'Units','Inches');
pos = get(fig3, 'Position');
set(fig3,'PaperPositionMode','Auto', 'PaperUnits','Inches', 'PaperSize',[pos(3),pos(4)]);
print(fig3, fullfile(cd,'Summary_LinearityAnalysis_BlandAltman_Position.pdf'), '-vector', '-dpdf');
savefig(fig3, 'Summary_LinearityAnalysis_BlandAltman_Position.fig');

%% Save data
saveAnalysisInfo_APP; 

end