function radialDifferences = computeRadialDifferences(image, Px, Py)
    %COMPUTERADIALDIFFERENCES Computes radial differences from (Px,Py).
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