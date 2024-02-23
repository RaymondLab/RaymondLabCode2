function newStack = preprocessImages(stack, enchanceContrast, roi)
    %PREPROCESSIMAGES
    
    % Upsample image for better detection of CR
    newStack = imresize(stack, 2);

    % Gaussian Filter
    newStack = imfilter(newStack, fspecial('gaussian',3,0.5));

    % Enhance Image Contrast
    % (cannot be done to stack, must loop)
    if enchanceContrast
        for i = 1:size(newStack,3)
            newStack(:,:,i) = imadjust(newStack(:,:,i), enchanceContrast, []);
        end
    end
    
    % Select elipse ROI
    % (cannot be done to stack, must loop)
    if exist('roi', 'var')
        for i = 1:size(newStack,3)
            tempImage = newStack(:,:,i);
                tempImage(~roi) = nan;
            newStack(:,:,i) = tempImage;
        end
    end
 
end