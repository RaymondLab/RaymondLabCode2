function [saccDilMask, saccMask, x_desaccaded] = run_saccadeDetection(x, base, thresh, lrPad, minChunkLen, method)

% Set default method if none is provided
if nargin<6, method = "SVT"; end

xmse = (x - base);
if strcmp(method, "SVT")
    saccMask = xmse.^2 > thresh;
else
    error('Invalid saccade detection "%s" was provided! Aborting script...', method);
end

% Create filter/kernel to dilate saccade regions in saccMask by "lrpad" samples on each side
padfilter = ones(1, 2*lrPad + 1);

% Convolution replaces zeros with ones around saccade regions (a.k.a. dilation)
dilatedMask = conv(double(saccMask), padfilter, 'same');

% x with desaccade segments removed
x_desaccaded = x;
x_desaccaded(logical(dilatedMask)) = NaN;

% Further refine saccDilMask to remove non-saccade segments that are "too small"
% Find where x_desaccaded is NaN (saccades)
[nanLocs_all, ~] = find(isnan(x_desaccaded));

% Find length of all 'good' (non-nan) chunks of data
goodChunk_len = nanLocs_all(2:end) - nanLocs_all(1:end-1);
goodChunk_len(goodChunk_len==1) = [];

% Find start times of all 'good' chunks
nanEnds_absolute = nanLocs_all(find((nanLocs_all(2:end)-nanLocs_all(1:end-1))-1));
goodChunk_starts = nanEnds_absolute + 1;

% Remove good chunks that are not too small
goodChunk_starts(goodChunk_len >= minChunkLen) = [];
goodChunk_len(goodChunk_len >= minChunkLen) = [];

% All good chunks that are smaller than specified length are NaN-ed
for kk = 1:length(goodChunk_starts)        
    x_desaccaded(goodChunk_starts(kk):goodChunk_starts(kk)+goodChunk_len(kk)) = NaN;
end

% Corresponding dilated and closed saccade mask
saccDilMask = isnan(x_desaccaded);

end