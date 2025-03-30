function [plots] = plotEyeTrackingImageFrame_APP(app, iImg, iFrameData, iPlotData, plots)
%PLOTEYETRACKINGIMAGEFRAME_APP Plots the i-th image and eye tracking result.

%% Initialize variables
theta = 0:0.1:2*pi;
phi   = linspace(0, 2*pi, 100);

%% Prepare x- and y-values for plots
cr1Lx = iFrameData.cr1_r.*cos(theta) + iFrameData.cr1_x;
cr1Ly = iFrameData.cr1_r.*sin(theta) + iFrameData.cr1_y;
cr1Mx = iFrameData.cr1_x;
cr1My = iFrameData.cr1_y;

cr2Lx = iFrameData.cr2_r.*cos(theta) + iFrameData.cr2_x;
cr2Ly = iFrameData.cr2_r.*sin(theta) + iFrameData.cr2_y;
cr2Mx = iFrameData.cr2_x;
cr2My = iFrameData.cr2_y;

pLx   = iFrameData.pupil_r1*cos(phi)*cos(iFrameData.pupil_angle) - ...
        sin(iFrameData.pupil_angle)*iFrameData.pupil_r2*sin(phi) + ...
        iFrameData.pupil_x;
pLy   = iFrameData.pupil_r1*cos(phi)*sin(iFrameData.pupil_angle) + ...
        cos(iFrameData.pupil_angle)*iFrameData.pupil_r2*sin(phi) + ...
        iFrameData.pupil_y;
pMx   = iFrameData.pupil_x;
pMy   = iFrameData.pupil_y;

if ~exist('plots', 'var')
    %% Generate new eye tracking plot
    cla(app.UIAxesTab1EyeTracking);
    plots.imgPlot = imagesc(app.UIAxesTab1EyeTracking, iImg);
    colormap(app.UIAxesTab1EyeTracking, 'gray');
    axis(app.UIAxesTab1EyeTracking, 'image');
    set(app.UIAxesTab1EyeTracking, 'ydir', 'reverse');
    set(app.UIAxesTab1EyeTracking, 'Visible', 'off');
    
    aaa = find(nanmean(iImg,1));
    bbb = find(nanmean(iImg,2));
    xlim(app.UIAxesTab1EyeTracking, [min(aaa),max(aaa)]);
    ylim(app.UIAxesTab1EyeTracking, [min(bbb),max(bbb)]);
    hold(app.UIAxesTab1EyeTracking, 'on');
    
    % Corneal Reflection 1
    plots.plot1 = line(app.UIAxesTab1EyeTracking, cr1Lx, cr1Ly, 'Color','b');
    plots.plot4 = plot(app.UIAxesTab1EyeTracking, cr1Mx, cr1My, '+b', 'LineWidth',2, 'MarkerSize',10);
    
    % Corneal Reflection 2
    plots.plot2 = line(app.UIAxesTab1EyeTracking, cr2Lx, cr2Ly, 'Color','c');
    plots.plot5 = plot(app.UIAxesTab1EyeTracking, cr2Mx, cr2My, '+c', 'LineWidth',2, 'MarkerSize',10);
    
    % Pupil
    plots.plot3 = line(app.UIAxesTab1EyeTracking, pLx, pLy, 'Color','r');
    plots.plot6 = plot(app.UIAxesTab1EyeTracking, pMx, pMy, '+r', 'LineWidth',2, 'MarkerSize',10);
    
    if exist('iPlotData','var')
        plots.plot7 = plot(app.UIAxesTab1EyeTracking, iPlotData.epx, iPlotData.epy, '.c');
        plots.plot8 = plot(app.UIAxesTab1EyeTracking, iPlotData.epx2, iPlotData.epy2, '.y');
    end
    hold(app.UIAxesTab1EyeTracking, 'off');
else
    %% Efficiently update existing eye tracking plot
    set(plots.imgPlot, 'CData',iImg);
    set(plots.plot1, 'XData',cr1Lx, 'YData',cr1Ly);
    set(plots.plot2, 'XData',cr2Lx, 'YData',cr2Ly);
    set(plots.plot3, 'XData',pLx, 'YData',pLy);
    set(plots.plot4, 'XData',cr1Mx, 'YData',cr1My);
    set(plots.plot5, 'XData',cr2Mx, 'YData',cr2My);
    set(plots.plot6, 'XData',pMx, 'YData',pMy);
    if exist('iPlotData','var')
        set(plots.plot7, 'XData',iPlotData.epx, 'YData',iPlotData.epy);
        set(plots.plot8, 'XData',iPlotData.epx2, 'YData',iPlotData.epy2);
    end
end

end