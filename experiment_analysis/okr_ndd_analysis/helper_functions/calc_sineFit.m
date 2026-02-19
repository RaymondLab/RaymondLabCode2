function [result, stat] = calc_sineFit(freq, samplerate, mask, eyevel, chairvel, drumvel)
%CALC_SINEFIT Calculates Least Squares Sinusoidal Fit of Provided Velocity Data
%   Detailed explanation goes here

% Validate necessary input parameters
if ~isnumeric(freq) || isempty(freq), error('Invalid "freq" was provided! Aborting calc_sineFit...'); end
if ~isnumeric(samplerate) || isempty(samplerate), error('Invalid "samplerate" was provided! Aborting calc_sineFit...'); end
if ~isnumeric(eyevel) || isempty(eyevel), error('Invalid "eyevel" was provided! Aborting calc_sineFit...'); end

% Initialize output struct
result = struct;

% Get time vector
nsamples = length(eyevel);
times    = (0:nsamples-1) / samplerate;

% Set up sinusoidal regressor matrix
regressors = [sin(2*pi*freq*times(:)), ...
              cos(2*pi*freq*times(:)), ...
              ones(nsamples,1)];

% Fit data using a least-squares linear regression
% b=coefficients, bint=95% confidence intervals, r=residual, rint=interval
% Stats include: 1) R-square, 2) F-statistics, 3) p-value 4) error variance
warning('off', 'stats:regress:NoConst') % No constant term bc sine centered at 0
cleanup = onCleanup(@() warning('on', 'stats:regress:NoConst'));

% ------------ EYE VELOCITY------------
eyevel = eyevel(:);
if ~isempty(mask)
    [b,~,~,~,stat] = regress(eyevel(mask), regressors(mask,:));
else
    [b,~,~,~,stat] = regress(eyevel, regressors);
end
eyevel_amp     = sqrt(b(1)^2 + b(2)^2);
eyevel_phase   = rad2deg(atan2(b(2), b(1)));
eyevel_offset  = b(3);
eyevel_fit     = regressors * b;

% Populate results
result.eyevel_amp    = eyevel_amp;
result.eyevel_phase  = eyevel_phase;
result.eyevel_offset = eyevel_offset;
result.eyevel_fit    = eyevel_fit;

% Regress chairvel and drumvel only if they are both provided
if nargin > 5
    % ------------CHAIR VELOCITY------------
    b               = regressors \ chairvel(:);  % Equivalent to: regress(chairvel(:), regressors);
    chairvel_amp    = sqrt(b(1)^2 + b(2)^2);
    chairvel_angle  = rad2deg(atan2(b(2), b(1)));
    chairvel_offset = b(3);
    
    % ------------DRUM VELOCITY------------
    b              = regressors \ drumvel(:);  % Equivalent to: regress(drumvel(:), regressors);
    drumvel_amp    = sqrt(b(1)^2 + b(2)^2);
    drumvel_angle  = rad2deg(atan2(b(2), b(1)));
    drumvel_offset = b(3);
    
    % Calculate eye gain as VOR/OKR based on chair/drum signal
    if chairvel_amp > 3     % Chair signal
        ref_amp   = chairvel_amp;
        ref_angle = chairvel_angle;        
    elseif drumvel_amp > 3  % Drum signal (no chair signal)
        ref_amp   = drumvel_amp;
        ref_angle = drumvel_angle;        
    else                    % No stimulus
        ref_amp   = 1;
        ref_angle = 0;
    end
    
    % ------------EYE RELATIVE TO CHAIR/DRUM------------
    eyevel_rel_gain  = eyevel_amp / ref_amp;
    eyevel_rel_phase = (eyevel_phase - ref_angle);
    eyevel_rel_phase = mod(eyevel_rel_phase, 360) - 180;
    eyevel_rel_fit   = sin(2*pi*freq*times + deg2rad(eyevel_rel_phase+180)) * eyevel_amp;

    % Add to results output
    result.chairvel_amp    = chairvel_amp;
    result.chairvel_angle  = chairvel_angle;
    result.chairvel_offset = chairvel_offset;
    result.drumvel_amp    = drumvel_amp;
    result.drumvel_angle  = drumvel_angle;
    result.drumvel_offset = drumvel_offset;
    result.ref_amp   = ref_amp;
    result.ref_angle = ref_angle;
    result.eyevel_rel_gain  = eyevel_rel_gain;
    result.eyevel_rel_phase = eyevel_rel_phase;
    result.eyevel_rel_fit   = eyevel_rel_fit(:);
end


end