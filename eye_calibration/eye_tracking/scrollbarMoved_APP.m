function scrollbarMoved_APP(source, event, app, trackParams, frameData)
%SCROLLBARMOVED

if length(event) > 1
    source.Position = event;
    currentPos = event(1,1);
else
    currentPos = round(event.CurrentPosition(1,1));
end
source.Label = sprintf(['Frame ', num2str(currentPos)]);
img = readImg(trackParams.defaultImgAdj, trackParams.cam, currentPos);
plotEyeTrackingImageFrame_APP(app, img, frameData(currentPos));
pause(0.01);  % Allow time for plot to update

end