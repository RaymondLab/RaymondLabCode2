function [cycleMat, cycleMean] = VOR_breakTrace(cycleLength, startpt, vector)

[cycleMat, ~] = vec2mat(vector(startpt:end), cycleLength, NaN);

% Remove final (partial) cycle
cycleMat(end,:) = [];

% Calculate Cycles Mean
cycleMean = nanmean(cycleMat, 1);