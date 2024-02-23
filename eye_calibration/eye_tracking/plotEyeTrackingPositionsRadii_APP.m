function [plots] = plotEyeTrackingPositionsRadii_APP(app, frameData, plots)
%PLOTEYETRACKINGPOSITIONSRADII_APP Plots eye tracking position and radii traces.

%% Initialize variables
timeLength = length([frameData.time1]);

%% Prepare y-axis limits for Radii plot
pupilRadii = max([frameData.pupil_r1], [frameData.pupil_r1]);
[minpr, maxpr] = bounds(pupilRadii);
[mincr1r, maxcr1r] = bounds([frameData.cr1_r]);
[mincr2r, maxcr2r] = bounds([frameData.cr2_r]);
minr = min([minpr, mincr1r, mincr2r])*0.85;
maxr = max([maxpr, maxcr1r, maxcr2r])*1.15;

if ~exist('plots', 'var')
    %% Generate new position and radii plots
    cla(app.UIAxesTab1PupilPos);
    cla(app.UIAxesTab1CRPos);
    cla(app.UIAxesTab1Radii);
    
    % Pupil position
    plots.plot1 = plot(app.UIAxesTab1PupilPos, ...
                       [frameData.pupil_x] - nanmean([frameData.pupil_x]), ...
                       'r-');
    ylabel(app.UIAxesTab1PupilPos, 'Pupil Position');
    yticks(app.UIAxesTab1PupilPos, -15:5:15);
    grid(app.UIAxesTab1PupilPos, 'on');
    set(app.UIAxesTab1PupilPos, 'GridAlpha',0.07);
    xlim(app.UIAxesTab1PupilPos, [0,timeLength]);
    ylim(app.UIAxesTab1PupilPos, [-20,20]);
    
    % CR1 and CR2 positions
    plots.plot2 = plot(app.UIAxesTab1CRPos, ...
                       [frameData.cr1_x] - nanmean([frameData.cr1_x]), ...
                       'b-');
    hold(app.UIAxesTab1CRPos, 'on');
    plots.plot3 = plot(app.UIAxesTab1CRPos, ...
                       [frameData.cr2_x] - nanmean([frameData.cr2_x]), ...
                       'c-');
    ylabel(app.UIAxesTab1CRPos, 'CR Position');
    yticks(app.UIAxesTab1CRPos, -15:5:15);
    grid(app.UIAxesTab1CRPos, 'on');
    set(app.UIAxesTab1CRPos, 'GridAlpha',0.07);
    xlim(app.UIAxesTab1CRPos, [0,timeLength]);
    ylim(app.UIAxesTab1CRPos, [-20,20]);
    hold(app.UIAxesTab1CRPos, 'off');
    
    % Pupil, CR1, and CR2 radii
    plots.plot4 = plot(app.UIAxesTab1Radii, pupilRadii, 'r-');
    hold(app.UIAxesTab1Radii, 'on');
    plots.plot5 = plot(app.UIAxesTab1Radii, [frameData.cr1_r],'b-');
    plots.plot6 = plot(app.UIAxesTab1Radii, [frameData.cr2_r],'c-');
    ylabel(app.UIAxesTab1Radii, 'Radii');
    xlim(app.UIAxesTab1Radii, [0,timeLength]);
    ylim(app.UIAxesTab1Radii, [minr,maxr]);
    hold(app.UIAxesTab1Radii, 'off');
    
    % Link the x-axes of Pupil and CR1/CR2 position plots
    linkaxes([app.UIAxesTab1PupilPos, app.UIAxesTab1CRPos, app.UIAxesTab1Radii], 'x');
else
    %% Efficiently update existing position and radii plots
    set(plots.plot1, 'YData', [frameData.pupil_x] - nanmean([frameData.pupil_x]));
    set(plots.plot2, 'YData', [frameData.cr1_x] - nanmean([frameData.cr1_x]));
    set(plots.plot3, 'YData', [frameData.cr2_x] - nanmean([frameData.cr2_x]));
    set(plots.plot4, 'YData', pupilRadii);
    set(plots.plot5, 'YData', [frameData.cr1_r]);
    set(plots.plot6, 'YData', [frameData.cr2_r]);
    ylim(app.UIAxesTab1Radii, [minr,maxr]);
end

end