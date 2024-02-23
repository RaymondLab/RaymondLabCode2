function saveplotMagnetVideoAlignment_APP(app)
%SAVEPLOTMAGNETVIDEOALIGNMENT_APP

%% Figure 1: Aligned Position and Velocity Traces
fig1 = figure('units','normalized', 'outerposition',[0 0 1 1]); clf
h1 = tiledlayout(fig1, 2, 1, 'TileSpacing','tight', 'Padding','compact');
ax1 = nexttile;
copyobj(app.UIAxesAlignedPositions.Children, ax1);
legend(ax1);
ax1.XLim = app.UIAxesAlignedPositions.XLim;
ax1.YLim = app.UIAxesAlignedPositions.YLim;
ax1.XLabel.String = app.UIAxesAlignedPositions.XLabel.String;
ax1.YLabel.String = app.UIAxesAlignedPositions.YLabel.String;
title(ax1, 'Alignment of Position Traces', 'FontSize',16);
grid(ax1, 'on');
box(ax1, 'on');

ax2 = nexttile;
copyobj(app.UIAxesAlignedVelocities.Children, ax2);
ax2.XLim = app.UIAxesAlignedVelocities.XLim;
ax2.YLim = app.UIAxesAlignedVelocities.YLim;
ax2.XLabel.String = app.UIAxesAlignedVelocities.XLabel.String;
ax2.YLabel.String = app.UIAxesAlignedVelocities.YLabel.String;
title(ax2, 'Alignment of Velocity Traces', 'FontSize',16);
grid(ax2, 'on');
box(ax2, 'on');

% Link the x-axes
linkaxes([ax1, ax2], 'x');

%% Figure 2: Aligned Velocity Cycles
fig2 = figure('units','normalized', 'outerposition',[0 0 1 1]);
axs2 = axes(fig2);
copyobj(app.UIAxesAlignedVelocityCycles.Children, axs2);
legend(axs2);
axs2.XLim = app.UIAxesAlignedVelocityCycles.XLim;
axs2.YLim = app.UIAxesAlignedVelocityCycles.YLim;
axs2.XLabel.String = app.UIAxesAlignedVelocityCycles.XLabel.String;
axs2.YLabel.String = app.UIAxesAlignedVelocityCycles.YLabel.String;
axs2.Title.String = app.UIAxesAlignedVelocityCycles.Title.String;
axs2.Title.FontSize = 16;
grid(axs2, 'on');
box(axs2, 'on');

%% Save Figures
set(fig1, 'Units','Inches');
pos = get(fig1, 'Position');
set(fig1,'PaperPositionMode','Auto', 'PaperUnits','Inches', 'PaperSize',[pos(3),pos(4)]);
print(fig1, fullfile(cd,'MagnetVideoAlignment_Traces.pdf'), '-vector', '-dpdf');
savefig('MagnetVideoAlignment_Traces.fig');

set(fig2, 'Units','Inches');
pos = get(fig2, 'Position');
set(fig2,'PaperPositionMode','Auto', 'PaperUnits','Inches', 'PaperSize',[pos(3),pos(4)]);
print(fig2, fullfile(cd,'MagnetVideoAlignment_Cycle.pdf'), '-vector', '-dpdf');
savefig('MagnetVideoAlignment_Cycle.fig');

close(fig1);
close(fig2);

end