% Camera and frame number of image to load
cam = 1;
frame = 1;

% XY-coordinate of pupil center in loaded image
xy = [288, 168];

% imageAdjust low/high parameters
lowIn = 0.1; %0.1;
highIn = 0.85; %0.85;

% Option to remove or not remove bright spots
removeBrightSpots = true;

% Percent of maximum intensity to remove
removeBrightSpotsPct = 0.8;

% Size of the disk for bright spot removal
removeBrightSpotsSize = 3;

% Colorbar limit
climit = 30;


%% Script starts here

% Load sample image
img = imread(fullfile(cd,['img',num2str(cam),'.tiff']),'Index',frame);

% Upsample image
img = imresize(img, 2);

% Contrast adjusted and normalized cubic transformed images
imgAdjust = imadjust(img,[lowIn,highIn],[]);
imgPupil = double(imgAdjust);
InoCR = imgPupil;

if removeBrightSpots
    [A, circen, crr] = CircularHough_Grd_TEST(imgAdjust, [6, 13], 10, 8, 1);
    maxvals = diag(imgAdjust(round(circen(:,2)),round(circen(:,1))));
    
    % Remove any bright spots
    removeCRthresh = max(maxvals)*removeBrightSpotsPct;
    disp(max(maxvals));
    totalMask = imgPupil > 204; %removeCRthresh;
    
    % NOTE: strel is what controls the size of the removed bright spot areas
    totalMaskDilate = imdilate(totalMask, strel('disk',removeBrightSpotsSize));
    InoCR(totalMaskDilate) = NaN;
end

% Apply gaussian smoothing to image
gaussian_smooth_image = @(I, sigma) imfilter(I, fspecial('gaussian', [ceil(2.5*sigma) ceil(2.5*sigma)], sigma), 'symmetric');
InoCR = gaussian_smooth_image(InoCR, 3);

% Compute Radial Differences of image
imgRadialDiff = computeRadialDifferences(InoCR, xy(1), xy(2));

dispMinMaxText('\nimgAdjust', imgAdjust);
dispMinMaxText('InoCR', InoCR);
dispMinMaxText('imgRadialDiff', imgRadialDiff);

%% Plot figure
fig = figure('units','normalized', 'outerposition',[0.05 0.2 0.9 0.6]);

ax(1) = subplot(1, 2, 1, 'Parent',fig);
ax1 = imagesc(ax(1), InoCR);
colormap(ax(1), 'gray');
set(ax(1), 'ydir', 'reverse');
set(ax(1), 'Visible', 'off');

ax(2) = subplot(1, 2, 2, 'Parent',fig);
ax2 = imagesc(ax(2), imgRadialDiff);
h2 = drawcrosshair(Position=xy, Color='w', LineWidth=0.5);
colormap(ax(2), 'jet');
colorbar(ax(2));
clim(ax(2), [0, climit]);
set(ax(2), 'ydir', 'reverse');
set(ax(2), 'Visible', 'off');

linkaxes([ax(1), ax(2)], 'xy');



function dispMinMaxText(name, image)
    minmaxtext = sprintf([name,' minmax = (',num2str(min(image(:))),', ', num2str(max(image(:))),')']);
    disp(minmaxtext);
end


function radialDifferences = computeRadialDifferences(image, Px, Py)
    % image: A 2D array of grayscale intensities
    % Px, Py: Coordinates of point P
    % Initialize the output array with zeros
    radialDifferences = zeros(size(image));
    
    % Get image dimensions
    [rows, cols] = size(image);
    
    % Iterate over each pixel in the image
    for x = 1:cols
        for y = 1:rows
            % Calculate the direction vector from P to the current pixel
            dx = x - Px;
            dy = y - Py;
            
            % Normalize the direction vector
            norm = sqrt(dx^2 + dy^2);
            if norm ~= 0
                dx = dx / norm;
                dy = dy / norm;
            end
            
            % Find the adjacent pixel coordinates in the radial direction towards P
            adjX = round(x - dx);
            adjY = round(y - dy);
            
            % Check if the adjacent pixel is within image bounds
            if adjX >= 1 && adjX <= cols && adjY >= 1 && adjY <= rows
                % Calculate the radial difference
                radialDifferences(y, x) = image(y, x) - image(adjY, adjX);
            end
        end
    end
end