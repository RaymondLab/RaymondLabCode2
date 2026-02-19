function cyclemat = segmentIntoCycles(data, startid, cyclelength)
%SEGMENTINTOCYCLES Updated version of legacy VOR_breakTrace function
%   Note: Per MATLAB recommendation, vec2mat was replaced with reshape.

% Segment data to start at first positive cycle
subdata = data(startid:end);  

% Corresponding number of sample columns and cycle rows
nsamples = numel(subdata);
ncycles = ceil(nsamples / cyclelength);  

% Pad with NaN
subdata(end+1:ncycles*cyclelength) = NaN;  
cyclemat = reshape(subdata, cyclelength, []).';

% Remove final cycle if is a partial cycle
npadded = ncycles * cyclelength - nsamples;
if npadded > 0
    cyclemat(end,:) = [];
end

end