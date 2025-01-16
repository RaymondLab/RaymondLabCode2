function vars = saveplotLinearityAnalysis_APP(app, vars)
%SAVEPLOTLINEARITYANALYSIS_APP

%% Load data
loadAnalysisInfo_APP;

%% Linearity Calculation
c = linspace(1, 10, length(mag1.pos_data_aligned_scaledInVel));

% Get saccade start and end times, exit early upon error
[mag1GoodStarts,mag1GoodStops,errInfo1] = getSaccadeStartEndTimes(mag1.saccades_all);
[mag2GoodStarts,mag2GoodStops,errInfo2] = getSaccadeStartEndTimes(mag1.saccades_all);
if errInfo1 || errInfo2
    return
end

% Get aligned desaccaded position data
mag1_aligned = mag1.pos_data_aligned_scaledInVel;
mag1_aligned(mag1.saccades_all) = nan;
vid1_aligned = vid.pos_data_upsampled_aligned;
vid1_aligned(mag1.saccades_all) = nan;

mag2_aligned = mag2.pos_data_aligned_scaledInVel;
mag2_aligned(mag2.saccades_all) = nan;
vid2_aligned = vid.pos_data_upsampled_aligned;
vid2_aligned(mag2.saccades_all) = nan;

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


%% Figure: Linearity of Each Magnet Channel
fig1 = figure('units','normalized', 'outerposition',[0 0 1 1]); clf
h1 = tiledlayout(fig1, 2, 3, 'TileSpacing','tight', 'Padding','compact');

%% Magnet Channel 1, Fit Coefficients
ax1 = nexttile;
for i = 1:length(mag1GoodStarts)
    chunk = mag1GoodStarts(i):mag1GoodStops(i);

    magPoints = mag1_aligned(chunk);
    magPoints = magPoints - nanmean(magPoints);
    mag1_aligned(chunk) = magPoints;
    
    vidChunk = vid1_aligned(chunk);
    vidChunk = vidChunk - nanmean(vidChunk);
    vid1_aligned(chunk) = vidChunk;
end

scatter(ax1, mag1_aligned, vid1_aligned, 50, c, '.', 'HandleVisibility','off');
hold(ax1, 'on');
xylimit = max(max(abs(mag1_aligned)), max(abs(vid1_aligned)));

slopes = [];
for i = 1:length(mag1GoodStarts)
    chunk = mag1GoodStarts(i)+1:mag1GoodStops(i);
    
    magPoints = mag1_aligned(chunk);
    magPoints = magPoints - nanmean(magPoints);
    
    vidChunk = vid1_aligned(chunk);
    vidChunk = vidChunk - nanmean(vidChunk);
    coefficients = polyfit(magPoints, vidChunk, 1);
    slopes(end+1) = coefficients(1);
    xFit = linspace(-xylimit, xylimit, 1000);
    yFit = polyval(coefficients , xFit);
    plot(ax1, xFit, yFit, '-k', 'LineWidth', .1, 'HandleVisibility','off');
end

coefficients = polyfit(mag1_aligned(~mag1.saccades_all), vid1_aligned(~mag1.saccades_all), 1);
xFit = linspace(-xylimit, xylimit, 1000);
yFit = polyval(coefficients , xFit);
plot(ax1, xFit, yFit, '-r', 'LineWidth',3, 'DisplayName',sprintf(' m = %.5f', coefficients(1)));
xline(ax1, 0, ':r', 'HandleVisibility','off');
yline(ax1, 0, ':r', 'HandleVisibility','off');
ylim(ax1, [-xylimit, xylimit]);
xlim(ax1, [-xylimit, xylimit]);
ax1.XLabel.String = 'Magnet Channel 1 Position';
ax1.YLabel.String = 'Video Channel Position';
title(ax1, 'Magnet Channel 1: Non-Saccade Segments', 'FontSize',14);
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
title(ax2, ['Magnet Channel 1 (with saccades): r^2: ', num2str(Rsq)], 'FontSize',14);
xlabel(ax2, sprintf('Video Position | minmax = (%.1f, %.1f)', vidmin, vidmax));
ylabel(ax2, sprintf('Magnet Channel 1 Position | minmax = (%.1f, %.1f)', magmin, magmax));
grid(ax2, 'on');
box(ax2, 'on');
hold(ax2, 'off');

ax3 = nexttile;
scatter(ax3, vid.pos_data_upsampled_aligned(~mag1.saccades_all), mag1.pos_data_aligned_scaledInVel(~mag1.saccades_all), 4, c(~mag1.saccades_all), '.');
[magmin, magmax] = bounds(mag1.pos_data_aligned_scaledInVel(~mag1.saccades_all));
[vidmin, vidmax] = bounds(vid.pos_data_upsampled_aligned(~mag1.saccades_all));
hold(ax3, 'on');
plot(ax3, mag1.pos_lin_all.yfit, mag1.pos_lin_all.range, 'k', 'lineWidth', 2);
colormap(ax3, hsv);
title(ax3, ['Magnet Channel 1 (saccades removed): r^2: ', num2str(mag1.pos_lin_all.Rsq)], 'FontSize',14)
xlabel(ax3, sprintf('Video Position | minmax = (%.1f, %.1f)', vidmin, vidmax));
ylabel(ax3, sprintf('Magnet Channel 1 Position | minmax = (%.1f, %.1f)', magmin, magmax));
grid(ax3, 'on');
box(ax3, 'on');
hold(ax3, 'off');
llimits = min(magmin, vidmin);
ulimits = max(magmax, vidmax);
ylim(ax3, [llimits, ulimits]);
xlim(ax3, [llimits, ulimits]);

%% Magnet Channel 2, Fit Coefficients
ax4 = nexttile;
for i = 1:length(mag2GoodStarts)
    chunk = mag2GoodStarts(i):mag2GoodStops(i);

    magPoints = mag2_aligned(chunk);
    magPoints = magPoints - nanmean(magPoints);
    mag2_aligned(chunk) = magPoints;
    
    vidChunk = vid2_aligned(chunk);
    vidChunk = vidChunk - nanmean(vidChunk);
    vid2_aligned(chunk) = vidChunk;
end

scatter(ax4, mag2_aligned, vid2_aligned, 50, c, '.', 'HandleVisibility','off');
hold(ax4, 'on');
xylimit = max(max(abs(mag2_aligned)), max(abs(vid2_aligned)));

for i = 1:length(mag2GoodStarts)
    chunk = mag2GoodStarts(i)+1:mag2GoodStops(i);
    
    magPoints = mag2_aligned(chunk);
    magPoints = magPoints - nanmean(magPoints);
    
    vidChunk = vid2_aligned(chunk);
    vidChunk = vidChunk - nanmean(vidChunk);
    coefficients = polyfit(magPoints, vidChunk, 1);
    xFit = linspace(-xylimit, xylimit, 1000);
    yFit = polyval(coefficients , xFit);
    plot(ax4, xFit, yFit, '-k', 'LineWidth', .1, 'HandleVisibility','off');
end

coefficients = polyfit(mag2_aligned(~mag2.saccades_all), vid2_aligned(~mag2.saccades_all), 1);
xFit = linspace(-xylimit, xylimit, 1000);
yFit = polyval(coefficients , xFit);
plot(ax4, xFit, yFit, '-r', 'LineWidth', 3, 'DisplayName',sprintf(' m = %.5f', coefficients(1)));
xline(ax4, 0, ':r', 'HandleVisibility','off');
yline(ax4, 0, ':r', 'HandleVisibility','off');
ylim(ax4, [-xylimit, xylimit]);
xlim(ax4, [-xylimit, xylimit]);
ax4.XLabel.String = 'Magnet Channel 2 Position';
ax4.YLabel.String = 'Video Channel Position';
title(ax4, 'Magnet Channel 2: Non-Saccade Segments', 'FontSize',14);
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
title(ax5, ['Magnet Channel 2 (with saccades): r^2: ', num2str(Rsq)], 'FontSize',14);
xlabel(ax5, sprintf('Video Position | minmax = (%.1f, %.1f)', vidmin, vidmax));
ylabel(ax5, sprintf('Magnet Channel 2 Position | minmax = (%.1f, %.1f)', magmin, magmax));
grid(ax5, 'on');
box(ax5, 'on');
hold(ax5, 'off');

ax6 = nexttile;
scatter(ax6, vid.pos_data_upsampled_aligned(~mag2.saccades_all), mag2.pos_data_aligned_scaledInVel(~mag2.saccades_all), 4, c(~mag2.saccades_all), '.');
[magmin, magmax] = bounds(mag2.pos_data_aligned_scaledInVel(~mag2.saccades_all));
[vidmin, vidmax] = bounds(vid.pos_data_upsampled_aligned(~mag2.saccades_all));
hold(ax6, 'on');
plot(ax6, mag2.pos_lin_all.yfit, mag2.pos_lin_all.range, 'k', 'lineWidth', 2);
colormap(ax6, hsv);
title(ax6, ['Magnet Channel 2 (saccades removed): r^2: ', num2str(mag2.pos_lin_all.Rsq)], 'FontSize',14)
xlabel(ax6, sprintf('Video Position | minmax = (%.1f, %.1f)', vidmin, vidmax));
ylabel(ax6, sprintf('Magnet Channel 2 Position | minmax = (%.1f, %.1f)', magmin, magmax));
grid(ax6, 'on');
box(ax6, 'on');
hold(ax6, 'off');
llimits = min(magmin, vidmin);
ulimits = max(magmax, vidmax);
ylim(ax6, [llimits, ulimits]);
xlim(ax6, [llimits, ulimits]);

%% Save Figures
set(fig1, 'Units','Inches');
pos = get(fig1, 'Position');
set(fig1,'PaperPositionMode','Auto', 'PaperUnits','Inches', 'PaperSize',[pos(3),pos(4)]);
print(fig1, fullfile(cd,'Summary_LinearityAnalysis.pdf'), '-vector', '-dpdf');
savefig('Summary_LinearityAnalysis.fig');

%% Save data
saveAnalysisInfo_APP;

end