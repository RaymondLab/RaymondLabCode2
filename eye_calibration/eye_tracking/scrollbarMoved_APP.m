function scrollbarMoved_APP(source, event, app, trackParams, frameData)
%SCROLLBARMOVED

if length(event) > 1
    source.Position = event;
    currentPos = event(1,1);
else
    currentPos = round(event.CurrentPosition(1,1));
end

lastIdx = round(trackParams.nImages - 1);
if currentPos < 1
    currentPos = 1;
    prevPos = event.PreviousPosition;
    event.CurrentPosition = [currentPos prevPos(2) currentPos prevPos(4)];
elseif currentPos > lastIdx
    currentPos = lastIdx;
    prevPos = event.PreviousPosition;
    event.CurrentPosition = [currentPos prevPos(2) currentPos prevPos(4)];
end

source.Label = sprintf(['Frame ', num2str(currentPos)]);
img = readImg(trackParams.defaultImgAdj, trackParams.cam, currentPos);
plotEyeTrackingImageFrame_APP(app, img, frameData(currentPos));
pause(0.01);  % Allow time for plot to update

end