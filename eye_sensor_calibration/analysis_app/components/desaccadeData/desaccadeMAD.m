function [saccadeLocs, eye_pos_filt, eye_vel_pfilt, threshMovMAD] = desaccadeMAD(data, samplerate, freq, presaccade, postsaccade, lambda, minDataLength, fc, mslen)
%DESACCADEMAD

% Apply lowpass Butterworth filter if needed
if ~exist('fc', 'var')
    % "data" assumed to be already filtered if "fc" is not provided
    eye_pos_filt = data;
else
    % "data" assumed to be raw position signal and is lowpass-filtered
    eye_pos_filt = butterworthfilter(data, fc, samplerate);
end
    
% Differentiate position to calculate velocity
veltau = .01;
if ~exist('mslen', 'var')
    if ~exist('fc', 'var')
        % "data" assumed to be filtered velocity if both "fc" and "mslen" are not provided
        eye_vel_pfilt = data;
    else
        eye_vel_pfilt = movingslope(eye_pos_filt, round(samplerate*veltau)) * samplerate;
    end
else
    eye_vel_pfilt = movingslope(eye_pos_filt, mslen) * samplerate;
end

% Remove experiment frequency by calculating the error from initial fit guess
segLength = length(data);
timeVec = 0:(1/samplerate):(segLength-1)/samplerate;

y1 = sin(2*pi*freq*timeVec(:));
y2 = cos(2*pi*freq*timeVec(:));
constant = ones(segLength,1);

vars = [y1, y2, constant];
keep = abs(eye_vel_pfilt) < 5*std(abs(eye_vel_pfilt)) + mean(abs(eye_vel_pfilt));
b = regress(eye_vel_pfilt(keep), vars(keep,:));
fit1 = vars *b;

eye_vel_err = (eye_vel_pfilt - fit1);

% Compute the moving MAD threshold
threshMovMAD = lambda*movmad(eye_vel_err, round(freq*samplerate));

% Remove all points that exceed the threshold
badDataLocations = abs(eye_vel_err - movmedian(eye_vel_err, round(freq*samplerate))) > threshMovMAD;

% Remove points around omit centers as defined by pre & post saccade time
presaccade = round((presaccade/1000)*samplerate);
postsaccade = round((postsaccade/1000)*samplerate);
sacmask = ones(1, presaccade+postsaccade);

% Filter function replaces zeros with ones (equal to remove time) around an omit center
rejecttemp1 = conv(double(badDataLocations), sacmask);
rejecttemp2 = rejecttemp1(presaccade:postsaccade+length(data)-1);

% Eyevel with desaccade segments removed
tempTrace = eye_pos_filt;
tempTrace(logical(rejecttemp2))= NaN;
saccadeLocs = isnan(tempTrace);

% Remove any 'good data' segments whose lengths are too small
goodStarts = strfind(saccadeLocs', [1 0]);
goodEnds = strfind(saccadeLocs', [0 1]);

if isempty(goodStarts)
    % BOTH are empty 
    if isempty(goodEnds)
        
    % START is empty
    else
        goodStarts = 1;
    end
else
    % END is empty
    if isempty(goodEnds)
        goodEnds = length(saccadeLocs);
        
    % NEITHER are empty
    else
        if goodEnds(1) < goodStarts(1)
            goodStarts = [1 goodStarts];
        end
        
        if goodEnds(end) < goodStarts(end)
            goodEnds = [goodEnds length(saccadeLocs)];
        end
    end
end  

for ii = 1:length(goodStarts)
    if (goodEnds(ii) - goodStarts(ii)) < minDataLength
        saccadeLocs(goodStarts(ii):goodEnds(ii)) = 1;
        goodStarts(ii) = NaN;
        goodEnds(ii) = NaN;
    end
end

end