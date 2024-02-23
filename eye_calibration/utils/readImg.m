function [img] = readImg(imadjVals, cam, frame, roi)
    %READIMG Single image processing pipeline for calibration.

    % Load sample images
    img = imread(fullfile(cd, ['img',num2str(cam),'.tiff']), 'Index',frame);
   
    % Upsample image for better detection of CR
    img = imresize(img, 2);

    % Enhance contrast
    img = imadjust(img, imadjVals, []);
    
    if exist('roi','var')
        img(~roi) = nan;
    end

end