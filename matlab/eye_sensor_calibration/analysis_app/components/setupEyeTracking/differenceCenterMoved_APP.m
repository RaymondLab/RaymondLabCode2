function differenceCenterMoved_APP(app, source, event)
%DIFFERENCECENTERMOVED_APP

img = app.data.rawImgStack(:,:,app.vars.frame);
updateTab2Plots(app, img);
pause(0.01);  % Allow time for plot to update

end