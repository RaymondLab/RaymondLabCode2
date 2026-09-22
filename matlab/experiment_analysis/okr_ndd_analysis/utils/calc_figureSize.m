function figure_size = calc_figureSize(resolution_scale)
%CALC_FIGURESIZE Calculates figure size based on resolution of display

% Set default resolution scale if not provided
if ~exist('resolution_scale', 'var')
    resolution_scale = 0.75;
end

% Calculate figure resolution based on size of smallest available display
monitors = get(groot, 'MonitorPositions');
[~, minIdx] = min(prod(monitors(:,3:4), 2));  % Index of smallest display
monxy = monitors(minIdx, 1:2);   % Origin (x,y) of the selected monitor
monwh = monitors(minIdx, 3:4);   % Width and height of the selected monitor
figwh = round(monwh * resolution_scale);
figxy = monxy + round((monwh - figwh) * 0.5);

% Calculated figure size
figure_size = [figxy, figwh];

end