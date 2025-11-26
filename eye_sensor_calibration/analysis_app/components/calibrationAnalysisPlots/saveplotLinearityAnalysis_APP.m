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
c1 = parula(length(mag1GoodStarts));
c2 = parula(length(mag2GoodStarts));

% Get aligned desaccaded position data
mag1_pos = mag1.pos_data_aligned_scaledInVel;
mag1_pos(mag1.saccades_all) = nan;
vid1_pos = vid.pos_data_upsampled_aligned;
vid1_pos(mag1.saccades_all) = nan;

mag2_pos = mag2.pos_data_aligned_scaledInVel;
mag2_pos(mag2.saccades_all) = nan;
vid2_pos = vid.pos_data_upsampled_aligned;
vid2_pos(mag2.saccades_all) = nan;

% Get aligned desaccaded velocity data
mag1_vel = mag1.vel_data_aligned_scaledInVel;
mag1_vel(mag1.saccades_all) = nan;
vid1_vel = vid.vel_data_upsampled_aligned;
vid1_vel(mag1.saccades_all) = nan;

mag2_vel = mag2.vel_data_aligned_scaledInVel;
mag2_vel(mag2.saccades_all) = nan;
vid2_vel = vid.vel_data_upsampled_aligned;
vid2_vel(mag2.saccades_all) = nan;

% Measure linearity of whole segments
mag1.pos_lin_all = measureLinearity(mag1.pos_data_aligned_scaledInVel(~mag1.saccades_all), vid.pos_data_upsampled_aligned(~mag1.saccades_all), ~mag1.saccades_all);
mag1.vel_lin_all = measureLinearity(mag1.vel_data_aligned_scaledInVel(~mag1.saccades_all), vid.vel_data_upsampled_aligned(~mag1.saccades_all), ~mag1.saccades_all);
mag2.pos_lin_all = measureLinearity(mag2.pos_data_aligned_scaledInVel(~mag2.saccades_all), vid.pos_data_upsampled_aligned(~mag2.saccades_all), ~mag2.saccades_all);
mag2.vel_lin_all = measureLinearity(mag2.vel_data_aligned_scaledInVel(~mag2.saccades_all), vid.vel_data_upsampled_aligned(~mag2.saccades_all), ~mag2.saccades_all);

% Measure linearity of non-saccades, split up into chunks
mag1_pos_zeromean = mag1_pos;
vid1_pos_zeromean = vid1_pos;
for i = 1:length(mag1GoodStarts)
    chunk = mag1GoodStarts(i):mag1GoodStops(i);
    magposchunk = mag1.pos_data_aligned_scaledInVel(chunk);
    vidposchunk = vid.pos_data_upsampled_aligned(chunk);
    magvelchunk = mag1.vel_data_aligned_scaledInVel(chunk);
    vidvelchunk = vid.vel_data_upsampled_aligned(chunk);
    mag1.pos_lin_chunks(i) = measureLinearity(magposchunk, vidposchunk);
    mag1.vel_lin_chunks(i) = measureLinearity(magvelchunk, vidvelchunk);
    mag1_pos_zeromean(chunk) = magposchunk - nanmean(magposchunk);
    vid1_pos_zeromean(chunk) = vidposchunk - nanmean(vidposchunk);
end

mag2_pos_zeromean = mag2_pos;
vid2_pos_zeromean = vid2_pos;
for i = 1:length(mag2GoodStarts)
    chunk = mag2GoodStarts(i):mag2GoodStops(i);
    magposchunk = mag2.pos_data_aligned_scaledInVel(chunk);
    vidposchunk = vid.pos_data_upsampled_aligned(chunk);
    magvelchunk = mag2.vel_data_aligned_scaledInVel(chunk);
    vidvelchunk = vid.vel_data_upsampled_aligned(chunk);
    mag2.pos_lin_chunks(i) = measureLinearity(magposchunk, vidposchunk);
    mag2.vel_lin_chunks(i) = measureLinearity(magvelchunk, vidvelchunk);
    mag2_pos_zeromean(chunk) = magposchunk - nanmean(magposchunk);
    vid2_pos_zeromean(chunk) = vidposchunk - nanmean(vidposchunk);
end

% Initialize data for Bland-Altman plots
xpavgs1 = 0.5 * (vid1_pos + mag1_pos);
ypdiffs1 = vid1_pos - mag1_pos;

xpavgs2 = 0.5 * (vid2_pos + mag2_pos);
ypdiffs2 = vid2_pos - mag2_pos;

p1 = getBlandAltmanParams(xpavgs1, ypdiffs1, mag1.saccades_all);
p2 = getBlandAltmanParams(xpavgs2, ypdiffs2, mag2.saccades_all);


%% Figure 1: Linearity of Position Traces of both Magnet Channels
fig1 = figure('units','normalized', 'outerposition',[0 0 1 1]); clf
h1 = tiledlayout(fig1, 2, 4, 'TileSpacing','tight', 'Padding','compact');

%% Magnet Channel 1 subplots
ax1 = nexttile;
ax2 = nexttile;
ax3 = nexttile;
ax8 = nexttile;

hold(ax1, 'on');
hold(ax2, 'on');
hold(ax3, 'on');
hold(ax8, 'on');

allData = logical(ones(length(~mag1.saccades_all), 1));

rsqs1 = [];
for i = 1:length(mag1GoodStarts)
    chunk = mag1GoodStarts(i)+1:mag1GoodStops(i);
    scatter(ax1, mag1_pos_zeromean(chunk), vid1_pos_zeromean(chunk), 50, c1(i,:), '.', 'HandleVisibility','off');
    scatter(ax3, mag1_pos(chunk), vid1_pos(chunk), 4, c1(i,:), '.', 'HandleVisibility','off');
    scatter(ax8, xpavgs1(chunk), ypdiffs1(chunk), 4, c1(i,:), '.', 'HandleVisibility','off');
    rsq_i = corrcoef(mag1_pos_zeromean(chunk), vid1_pos_zeromean(chunk));
    rsqs1(end+1) = rsq_i(1, 2).^2;
end
rsq1_mean = mean(rsqs1);
rsq1_med = median(rsqs1);

slopes = [];
xylimit = max(max(abs(mag1_pos_zeromean)), max(abs(vid1_pos_zeromean)));
coefficients = polyfit(mag1_pos_zeromean(~mag1.saccades_all), vid1_pos_zeromean(~mag1.saccades_all), 1);
for i = 1:length(mag1GoodStarts)
    chunk = mag1GoodStarts(i)+1:mag1GoodStops(i);
    coeffs = polyfit(mag1_pos_zeromean(chunk), vid1_pos_zeromean(chunk), 1);
    slopes(end+1) = coeffs(1) - coefficients(1);
    xFit = linspace(-xylimit, xylimit, 1000);
    yFit = polyval(coeffs , xFit);
    plot(ax1, xFit, yFit, '-k', 'LineWidth', .1, 'HandleVisibility','off');
end
mstd = std(slopes);
mvar = var(slopes);
xFit = linspace(-xylimit, xylimit, 1000);
yFit = polyval(coefficients , xFit);
plot(ax1, xFit, yFit, '-r', 'LineWidth',3, 'DisplayName',sprintf(' m = %.5f (std=%.3f, var=%.3f) | r^2 = %.5f (%.5f median)', coefficients(1), mstd, mvar, rsq1_mean, rsq1_med));
xline(ax1, 0, ':r', 'HandleVisibility','off');
yline(ax1, 0, ':r', 'HandleVisibility','off');
ylim(ax1, [-xylimit, xylimit]);
xlim(ax1, [-xylimit, xylimit]);
ax1.XLabel.String = 'Magnet Channel 1 Position';
ax1.YLabel.String = 'Video Channel Position';
title(ax1, sprintf('Magnet Channel 1 (saccades removed): (%d) Zero-Mean Segments',length(slopes)), 'FontSize',12);
grid(ax1, 'on');
box(ax1, 'on');
hold(ax1, 'off');
legend(ax1, 'Location','northwest', 'FontSize',11);

Rsq = linearityScatterPlot(ax2, mag1.pos_data_aligned_scaledInVel, vid.pos_data_upsampled_aligned, allData, c);
[magmin, magmax] = bounds(mag1.pos_data_aligned_scaledInVel(allData));
[vidmin, vidmax] = bounds(vid.pos_data_upsampled_aligned(allData));
title(ax2, sprintf('Magnet Channel 1 (with saccades): r^2 = %.5f', Rsq), 'FontSize',12);
xlabel(ax2, sprintf('Magnet Channel 1 Position | minmax = (%.1f, %.1f)', magmin, magmax));
ylabel(ax2, sprintf('Video Position | minmax = (%.1f, %.1f)', vidmin, vidmax));
grid(ax2, 'on');
box(ax2, 'on');
hold(ax2, 'off');

[magmin, magmax] = bounds(mag1_pos(~mag1.saccades_all));
[vidmin, vidmax] = bounds(vid1_pos(~mag1.saccades_all));
llimits = min(magmin, vidmin);
ulimits = max(magmax, vidmax);
plot(ax3, mag1.pos_lin_all.range, mag1.pos_lin_all.yfit, 'k', 'lineWidth', 2);
colormap(ax3, hsv);
title(ax3, sprintf('Magnet Channel 1 (saccades removed): r^2 = %.5f', mag1.pos_lin_all.Rsq), 'FontSize',12)
xlabel(ax3, sprintf('Magnet Channel 1 Position | minmax = (%.1f, %.1f)', magmin, magmax));
ylabel(ax3, sprintf('Video Position | minmax = (%.1f, %.1f)', vidmin, vidmax));
grid(ax3, 'on');
box(ax3, 'on');
ylim(ax3, [llimits, ulimits]);
xlim(ax3, [llimits, ulimits]);
hold(ax3, 'off');

yline(ax8, p1.b, '--k', sprintf('bias = %.4f',p1.b), 'LineWidth',1.25, 'FontSize',12, 'HandleVisibility','off');
yline(ax8, p1.UB, '--k', sprintf('UB = %.4f',p1.UB), 'LineWidth',1.25, 'FontSize',12, 'HandleVisibility','off');
yline(ax8, p1.LB, '--k', sprintf('LB = %.4f',p1.LB), 'LineWidth',1.25, 'FontSize',12, 'HandleVisibility','off');
plot(ax8, p1.xfit, p1.yfit1, '-k', 'LineWidth',3, 'DisplayName',p1.txt1);
plot(ax8, p1.xfit, p1.yfit3, '-r', 'LineWidth',2.5, 'DisplayName',p1.txt3);
xlabel(ax8, 'Mean of VIDEO and MAG1 positions ');
ylabel(ax8, 'VIDEO - MAG1 (position)');
xlim(ax8, p1.xlims);
ylim(ax8, p1.ylims);
hold(ax8, 'off');
grid(ax8, 'on');
box(ax8, 'on');
title(ax8, 'Magnet Channel 1 (saccades removed)', 'FontSize',12);
legend(ax8, 'Location','northwest', 'FontSize',10);


%% Magnet Channel 2 subplots
ax4 = nexttile;
ax5 = nexttile;
ax6 = nexttile;
ax10 = nexttile;

hold(ax4, 'on');
hold(ax5, 'on');
hold(ax6, 'on');
hold(ax10, 'on');

allData = logical(ones(length(~mag2.saccades_all), 1));

rsqs2 = [];
for i = 1:length(mag2GoodStarts)
    chunk = mag2GoodStarts(i)+1:mag2GoodStops(i);
    scatter(ax4, mag2_pos_zeromean(chunk), vid2_pos_zeromean(chunk), 50, c2(i,:), '.', 'HandleVisibility','off');
    scatter(ax6, mag2_pos(chunk), vid2_pos(chunk), 4, c2(i,:), '.', 'HandleVisibility','off');
    scatter(ax10, xpavgs2(chunk), ypdiffs2(chunk), 4, c2(i,:), '.', 'HandleVisibility','off');
    rsq_i = corrcoef(mag2_pos_zeromean(chunk), vid2_pos_zeromean(chunk));
    rsqs2(end+1) = rsq_i(1, 2).^2;
end
rsq2_mean = mean(rsqs2);
rsq2_med = median(rsqs2);

slopes = [];
xylimit = max(max(abs(mag2_pos_zeromean)), max(abs(vid2_pos_zeromean)));
coefficients = polyfit(mag2_pos_zeromean(~mag2.saccades_all), vid2_pos_zeromean(~mag2.saccades_all), 1);
for i = 1:length(mag2GoodStarts)
    chunk = mag2GoodStarts(i)+1:mag2GoodStops(i);
    coeffs = polyfit(mag2_pos_zeromean(chunk), vid2_pos_zeromean(chunk), 1);
    slopes(end+1) = coeffs(1) - coefficients(1);
    xFit = linspace(-xylimit, xylimit, 1000);
    yFit = polyval(coeffs , xFit);
    plot(ax4, xFit, yFit, '-k', 'LineWidth', .1, 'HandleVisibility','off');
end
mstd = nanstd(slopes);
mvar = nanvar(slopes);
xFit = linspace(-xylimit, xylimit, 1000);
yFit = polyval(coefficients , xFit);
plot(ax4, xFit, yFit, '-r', 'LineWidth', 3, 'DisplayName',sprintf(' m = %.5f (std=%.3f, var=%.3f) | r^2 = %.5f (%.5f median)', coefficients(1), mstd, mvar, rsq2_mean, rsq2_med));
xline(ax4, 0, ':r', 'HandleVisibility','off');
yline(ax4, 0, ':r', 'HandleVisibility','off');
ylim(ax4, [-xylimit, xylimit]);
xlim(ax4, [-xylimit, xylimit]);
ax4.XLabel.String = 'Magnet Channel 2 Position';
ax4.YLabel.String = 'Video Channel Position';
title(ax4, sprintf('Magnet Channel 2 (saccades removed): (%d) Zero-Mean Segments',length(slopes)), 'FontSize',12);
grid(ax4, 'on');
box(ax4, 'on');
hold(ax4, 'off');
legend(ax4, 'Location','northwest', 'FontSize',11);

Rsq = linearityScatterPlot(ax5, mag2.pos_data_aligned_scaledInVel, vid.pos_data_upsampled_aligned, allData, c);
[magmin, magmax] = bounds(mag2.pos_data_aligned_scaledInVel(allData));
[vidmin, vidmax] = bounds(vid.pos_data_upsampled_aligned(allData));
title(ax5, sprintf('Magnet Channel 2 (with saccades): r^2 = %.5f', Rsq), 'FontSize',12);
xlabel(ax5, sprintf('Magnet Channel 2 Position | minmax = (%.1f, %.1f)', magmin, magmax));
ylabel(ax5, sprintf('Video Position | minmax = (%.1f, %.1f)', vidmin, vidmax));
grid(ax5, 'on');
box(ax5, 'on');
hold(ax5, 'off');

[magmin, magmax] = bounds(mag2_pos(~mag2.saccades_all));
[vidmin, vidmax] = bounds(vid2_pos(~mag2.saccades_all));
llimits = min(magmin, vidmin);
ulimits = max(magmax, vidmax);
plot(ax6, mag2.pos_lin_all.range, mag2.pos_lin_all.yfit, 'k', 'lineWidth', 2);
colormap(ax6, hsv);
title(ax6, sprintf('Magnet Channel 2 (saccades removed): r^2 = %.5f', mag2.pos_lin_all.Rsq), 'FontSize',12)
xlabel(ax6, sprintf('Magnet Channel 2 Position | minmax = (%.1f, %.1f)', magmin, magmax));
ylabel(ax6, sprintf('Video Position | minmax = (%.1f, %.1f)', vidmin, vidmax));
grid(ax6, 'on');
box(ax6, 'on');
ylim(ax6, [llimits, ulimits]);
xlim(ax6, [llimits, ulimits]);
hold(ax6, 'off');

yline(ax10, p2.b, '--k', sprintf('bias = %.4f',p2.b), 'LineWidth',1.25, 'FontSize',12, 'HandleVisibility','off');
yline(ax10, p2.UB, '--k', sprintf('UB = %.4f',p2.UB), 'LineWidth',1.25, 'FontSize',12, 'HandleVisibility','off');
yline(ax10, p2.LB, '--k', sprintf('LB = %.4f',p2.LB), 'LineWidth',1.25, 'FontSize',12, 'HandleVisibility','off');
plot(ax10, p2.xfit, p2.yfit1, '-k', 'LineWidth',3, 'DisplayName',p2.txt1);
plot(ax10, p2.xfit, p2.yfit3, '-r', 'LineWidth',2.5, 'DisplayName',p2.txt3);
xlabel(ax10, 'Mean of VIDEO and MAG2 positions ');
ylabel(ax10, 'VIDEO - MAG2 (position)');
xlim(ax10, p2.xlims);
ylim(ax10, p2.ylims);
hold(ax10, 'off');
grid(ax10, 'on');
box(ax10, 'on');
title(ax10, 'Magnet Channel 2 (saccades removed)', 'FontSize',12);
legend(ax10, 'Location','northwest', 'FontSize',10);

linkaxes([ax1,ax4], 'xy');
linkaxes([ax2,ax5], 'xy');
linkaxes([ax3,ax6], 'xy');
linkaxes([ax8,ax10], 'xy');

title(h1, 'Linearity (Rainbow) Plots | Best Channel Highlighted with Blue Border', 'FontSize',18);

%% Add border around channel subplot with highest r^2 value
if rsq1_mean > rsq2_mean
    axpos = ax1.Position;
    vars.vid.bestMagChannel = 1;
else
    axpos = ax4.Position;
    vars.vid.bestMagChannel = 2;
end
annotation('rectangle', axpos, 'Color','blue', 'LineWidth',3);


%% Save Figure
set(fig1, 'Units','Inches');
pos = get(fig1, 'Position');
set(fig1,'PaperPositionMode','Auto', 'PaperUnits','Inches', 'PaperSize',[pos(3),pos(4)]);
print(fig1, fullfile(cd,'Summary_LinearityAnalysis_Position.pdf'), '-vector', '-dpdf');
savefig(fig1, 'Summary_LinearityAnalysis_Position.fig');

%% Save data
saveAnalysisInfo_APP; 

end