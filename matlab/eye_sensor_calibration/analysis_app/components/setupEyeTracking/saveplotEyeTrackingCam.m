function saveplotEyeTrackingCam(frameData, cam)
%SAVEPLOTEYETRACKINGCAM

fig = figure('units','normalized', 'outerposition',[0 0 1 1]); clf
h = tiledlayout(fig, 2, 1, 'TileSpacing','tight', 'Padding','compact');

ax1 = nexttile;
plot(ax1, [frameData.time1], [frameData.pupil_x], 'r-'); 
hold(ax1, 'on');
plot(ax1, [frameData.time1], [frameData.cr1_x], 'b-');
plot(ax1, [frameData.time2], [frameData.cr2_x], 'c-');
xticks(ax1, []);
xlim(ax1, [frameData(1).time1, frameData(end).time1]);
ylabel(ax1, 'Horizontal Position (px)');
grid(ax1, 'on');
box(ax1, 'on');
hold(ax1, 'off');

ax2 = nexttile;
plot(ax2, [frameData.time1], nanmax([frameData.pupil_r1; frameData.pupil_r2],[],1), 'r-');
hold(ax2, 'on');
plot([frameData.time1], [frameData.cr1_r], 'b-');
plot([frameData.time2], [frameData.cr2_r], 'c-');  
xlim(ax2, [frameData(1).time1, frameData(end).time1]);
xlabel(ax2, 'Time (s)');
ylabel(ax2, 'Radii (px)');
legend(ax2, {'Pupil', 'CR1', 'CR2'});
grid(ax2, 'on');
box(ax2, 'on');
hold(ax2, 'off');

set(fig, 'Units','Inches');
pos = get(fig, 'Position');
set(fig, 'PaperPositionMode','Auto', 'PaperUnits','Inches', 'PaperSize',[pos(3),pos(4)]);
print(fig, fullfile(cd,['Summary_EyeTracking_Cam',num2str(cam),'.pdf']), '-vector', '-dpdf');
savefig(['Summary_EyeTracking_Cam', num2str(cam), '.fig']);

pause(0.1);
close(fig);

end