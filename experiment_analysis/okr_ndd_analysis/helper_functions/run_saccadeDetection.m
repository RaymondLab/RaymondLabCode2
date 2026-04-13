function [saccDilMask, saccMask, x_desaccaded] = run_saccadeDetection(x, base, thresh, lrPad, minChunkLen, method)

% Set default method if none is provided
if nargin<6, method = "SVT"; end

xmse = (x - base);
if strcmp(method, "SVT")
    saccMask = xmse.^2 > thresh;
elseif strcmp(method, "MAD")
    maxIter     = 10;
    lambdaPeak  = thresh;  % Strict factor for onset/offset walk
    lambdaOnset = 3;       % Lenient factor for onset/offset walk
    % Start by assuming all finite samples are "clean"
    keep = ~isnan(xmse);  
    % Iteratatively re-compute threshold until it stabilizes
    for iter = 1:maxIter
        med                  = median(xmse(keep), 'omitnan');
	    sigma                = 1.4826 * median(abs(xmse(keep) - med), 'omitnan');
	    newKeep              = abs(xmse - med) <= thresh * sigma;
	    newKeep(isnan(xmse)) = false;  % NaNs never count as "clean"
	    % Stop iterating once the new threshold results in no new detections
	    if isequal(newKeep, keep), break, end
        keep = newKeep;
    end

    % Dual-threshold onset/offset walk
    peakMask  = abs(xmse - med) >  lambdaPeak  * sigma;   % strict
    onsetMask = abs(xmse - med) >  lambdaOnset * sigma;   % lenient
    onsetMask(isnan(xmse)) = false;

    % Keep only onsetMask regions that contain at least one peakMask sample
    CC = bwconncomp(onsetMask);
    saccMask = false(size(xmse));
    for k = 1:CC.NumObjects
        idx = CC.PixelIdxList{k};
        if any(peakMask(idx))
            saccMask(idx) = true;
        end
    end
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