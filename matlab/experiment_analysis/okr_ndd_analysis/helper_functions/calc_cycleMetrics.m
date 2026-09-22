function results = calc_cycleMetrics(hevel, stimvel)
% CALC_CYCLEMETRICS Computes VOR response characterizations per half-cycle
%
% SYNTAX:
%   results = calc_cycleMetrics(hevel, stimvel)
%
% INPUTS:
%   hevel   - Horizontal eye velocity array of size CYCLELENGTHxNCYCLES or 
%             NCYCLESxCYCLELENGTH (cycle length must be the larger dimension)
%   stimvel - Stimulus velocity array, must have the same dimensions as hevel
%
% OUTPUTS:
%   results - Struct containing characterization metrics organized by group.
%             Each metric is an NCYCLES x 2 matrix where:
%               Column 1 = first half-cycle
%               Column 2 = second half-cycle
%
% =========================================================================
% METRIC GROUP 1: CENTROID (MEAN) BASED
% =========================================================================
% These metrics use the temporal center of mass, which is analogous to the
% mean of a probability distribution. Sensitive to extreme values and tails.
%
%   Metric 1A - Temporal Centroid (Mean):
%       .eye.centroidMean   - Eye velocity center of mass time (ms)
%       .stim.centroidMean  - Stimulus velocity center of mass time (ms)
%
%       Interpretation: The "balance point" of the waveform in time. If the
%       waveform were a physical mass distribution along the time axis, this
%       is where it would balance. Earlier values indicate mass concentrated
%       toward the start of the half-cycle; later values indicate mass
%       concentrated toward the end.
%
%   Metric 1B - Lead/Lag via Centroid (Mean) Difference:
%       .centroidMeanDiff   - Eye centroidMean minus Stim centroidMean (ms)
%
%       Interpretation: The temporal shift between response and stimulus
%       centers of mass. Positive values indicate the eye response lags
%       the stimulus; negative values indicate the eye response leads.
%
%   Metric 1C - Temporal Skew from Centroid (Mean):
%       .eye.skewFromCentroidMean   - Eye velocity skew index
%       .stim.skewFromCentroidMean  - Stimulus velocity skew index
%
%       Formula: 1 - (2 * centroidMean / halfCycleDuration)
%
%       Interpretation: Measures where the center of mass falls relative to
%       the temporal midpoint of the half-cycle.
%         Positive: Mass concentrated early (waveform "leans right")
%         Zero:     Mass centered at midpoint (symmetric)
%         Negative: Mass concentrated late (waveform "leans left")
%
% =========================================================================
% METRIC GROUP 2: CENTROID (MEDIAN) AND QUANTILE BASED
% =========================================================================
% These metrics use cumulative area to find quantile times, analogous to
% percentiles of a probability distribution. Robust to extreme values.
%
%   Metric 2A - Temporal Centroid (Median):
%       .eye.centroidMedian   - Eye velocity 50th percentile time (ms)
%       .stim.centroidMedian  - Stimulus velocity 50th percentile time (ms)
%
%       Interpretation: The time at which 50% of the total response area
%       has occurred. Unlike the mean-based centroid, this is robust to
%       outlier peaks or long tails. Represents the "middle" of the response.
%
%   Metric 2B - Lead/Lag via Centroid (Median) Difference:
%       .centroidMedianDiff   - Eye centroidMedian minus Stim centroidMedian (ms)
%
%       Interpretation: Same as 1B but using the median-based centroid.
%       More robust to waveform distortions. Positive = lag, Negative = lead.
%
%   Metric 2C - Temporal Skew from Centroid (Median):
%       .eye.skewFromCentroidMedian   - Eye velocity skew index
%       .stim.skewFromCentroidMedian  - Stimulus velocity skew index
%
%       Formula: 1 - (2 * centroidMedian / halfCycleDuration)
%
%       Interpretation: Same as 1C but using the median-based centroid.
%       Positive = median early, Negative = median late.
%
%   Metric 2D - Quantile Times:
%       .eye.t25, .eye.t75    - Eye velocity 25th and 75th percentile times (ms)
%       .stim.t25, .stim.t75  - Stimulus velocity 25th and 75th percentile times (ms)
%
%       Interpretation: Time points at which 25% and 75% of total response
%       area has occurred. The difference (t75 - t25) is the interquartile
%       range (IQR), representing the temporal "spread" of the response.
%
%   Metric 2E - Quantile Skew (Bowley Skewness):
%       .eye.skewQuantile   - Eye velocity Bowley skewness
%       .stim.skewQuantile  - Stimulus velocity Bowley skewness
%
%       Formula: [(t75 - t50) - (t50 - t25)] / (t75 - t25)
%
%       Interpretation: Measures asymmetry in the spread around the median.
%         Positive: Upper half (t50 to t75) is wider than lower half
%                   (distribution has a late/right tail)
%         Zero:     Symmetric spread around median
%         Negative: Lower half (t25 to t50) is wider than upper half
%                   (distribution has an early/left tail)
%
% =========================================================================
% METRIC GROUP 3: MEAN-MEDIAN COMPARISON
% =========================================================================
% Comparing centroid (mean) to centroid (median) reveals the influence of
% tails and outliers on the waveform's central tendency.
%
%   Metric 3A - Mean-Median Skew:
%       .eye.skewMeanMedian   - Eye velocity mean-median comparison
%       .stim.skewMeanMedian  - Stimulus velocity mean-median comparison
%
%       Formula: (centroidMean - centroidMedian) / IQR
%       Where IQR = t75 - t25
%
%       Interpretation: Quantifies how much tails pull the mean away from
%       the median, normalized by the spread of the distribution.
%         Positive: Mean is later than median (late tail pulls mean rightward)
%         Zero:     Mean equals median (symmetric distribution)
%         Negative: Mean is earlier than median (early tail pulls mean leftward)
%
%       When eye and stim values differ substantially, it indicates the
%       response has different tail behavior than the stimulus.
%
% =========================================================================
% METRIC GROUP 4: PEAK-BASED
% =========================================================================
% These metrics identify the two tallest peaks in each half-cycle using
% findpeaks on the non-rectified eye velocity signal.
%
%   Metric 4A - Two Tallest Peaks:
%       .eye.peak1Val     - Amplitude of the tallest peak (deg/s)
%       .eye.peak1TimeMs  - Time of the tallest peak (ms)
%       .eye.peak2Val     - Amplitude of the second tallest peak (deg/s)
%       .eye.peak2TimeMs  - Time of the second tallest peak (ms)
%
%       Interpretation: A naive characterization of the dominant peaks in
%       the half-cycle response. If only one peak is detected, the second
%       peak fields are NaN. If no peaks are detected, both are NaN.
%       Peaks are found on the non-rectified signal with a minimum height of
%       0.01 deg/s and a minimum inter-peak distance of 50 ms.
%
% =========================================================================
%
% =========================================================================
%
% EXAMPLE:
%   % Generate sample data: 1000 samples/cycle, 10 cycles
%   t = (0:999)' / 1000;
%   stimvel = repmat(sin(2*pi*t), 1, 10);
%   hevel = repmat(sin(2*pi*t - 0.1), 1, 10) + 0.05*randn(1000, 10);
%   results = calc_cycleMetrics(hevel, stimvel);
%
%   % Compare lead/lag using mean vs median centroid
%   lagMean = results.centroidMeanDiff;
%   lagMedian = results.centroidMedianDiff;
%
% =========================================================================

    %% Input validation and orientation
    if nargin < 2
        error('Both hevel and stimvel are required inputs.');
    end
    
    [hevel, nCycles, cycleLength] = orientData(hevel);
    [stimvel, nCyclesStim, cycleLengthStim] = orientData(stimvel);
    
    % Validate dimensions match
    if nCycles ~= nCyclesStim || cycleLength ~= cycleLengthStim
        error('hevel and stimvel must have the same dimensions.');
    end
    
    %% Setup
    halfLength = cycleLength / 2;
    if mod(cycleLength, 2) ~= 0
        error('Cycle length must be even to split into equal half-cycles.');
    end
    
    % Time vector for half-cycle (in ms, assuming 1000 Hz sample rate)
    tHalf = (0:halfLength-1)';
    halfDuration = halfLength;
    
    %% Preallocate output arrays (NCYCLES x 2)
    
    % Group 1: Centroid (Mean) based
    results.eye.centroidMean = zeros(nCycles, 2);
    results.stim.centroidMean = zeros(nCycles, 2);
    results.centroidMeanDiff = zeros(nCycles, 2);
    results.eye.skewFromCentroidMean = zeros(nCycles, 2);
    results.stim.skewFromCentroidMean = zeros(nCycles, 2);
    
    % Group 2: Centroid (Median) and Quantile based
    results.eye.centroidMedian = zeros(nCycles, 2);
    results.stim.centroidMedian = zeros(nCycles, 2);
    results.centroidMedianDiff = zeros(nCycles, 2);
    results.eye.skewFromCentroidMedian = zeros(nCycles, 2);
    results.stim.skewFromCentroidMedian = zeros(nCycles, 2);
    results.eye.t25 = zeros(nCycles, 2);
    results.eye.t75 = zeros(nCycles, 2);
    results.stim.t25 = zeros(nCycles, 2);
    results.stim.t75 = zeros(nCycles, 2);
    results.eye.skewQuantile = zeros(nCycles, 2);
    results.stim.skewQuantile = zeros(nCycles, 2);
    
    % Group 3: Mean-Median comparison
    results.eye.skewMeanMedian = zeros(nCycles, 2);
    results.stim.skewMeanMedian = zeros(nCycles, 2);

    % Group 4: Peak-based
    results.eye.peak1Val = NaN(nCycles, 2);
    results.eye.peak1TimeMs = NaN(nCycles, 2);
    results.eye.peak2Val = NaN(nCycles, 2);
    results.eye.peak2TimeMs = NaN(nCycles, 2);
    
    %% Process each cycle
    for c = 1:nCycles
        % Extract current cycle
        eyeCycle = hevel(:, c);
        stimCycle = stimvel(:, c);
        
        % Split into half-cycles
        eyeHalves = {eyeCycle(1:halfLength), eyeCycle(halfLength+1:end)};
        stimHalves = {stimCycle(1:halfLength), stimCycle(halfLength+1:end)};
        
        for h = 1:2
            eyeHalf = eyeHalves{h};
            stimHalf = stimHalves{h};
            
            % Rectify signals for consistent computation across half-cycles
            eyeRect = abs(eyeHalf);
            stimRect = abs(stimHalf);
            
            %% Metric 1A: Temporal Centroid (Mean)
            eyeCentroidMean = computeCentroid(eyeRect, tHalf);
            stimCentroidMean = computeCentroid(stimRect, tHalf);
            results.eye.centroidMean(c, h) = eyeCentroidMean;
            results.stim.centroidMean(c, h) = stimCentroidMean;
            
            %% Metric 1B: Lead/Lag via Centroid (Mean) Difference
            results.centroidMeanDiff(c, h) = eyeCentroidMean - stimCentroidMean;
            
            %% Metric 1C: Temporal Skew from Centroid (Mean)
            results.eye.skewFromCentroidMean(c, h) = 1 - (2 * eyeCentroidMean / halfDuration);
            results.stim.skewFromCentroidMean(c, h) = 1 - (2 * stimCentroidMean / halfDuration);
            
            %% Metric 2A & 2D: Quantile Times (including Centroid Median = t50)
            [eyeT25, eyeT50, eyeT75] = computeQuantileTimes(eyeRect, tHalf);
            [stimT25, stimT50, stimT75] = computeQuantileTimes(stimRect, tHalf);
            
            results.eye.centroidMedian(c, h) = eyeT50;
            results.stim.centroidMedian(c, h) = stimT50;
            results.eye.t25(c, h) = eyeT25;
            results.eye.t75(c, h) = eyeT75;
            results.stim.t25(c, h) = stimT25;
            results.stim.t75(c, h) = stimT75;
            
            %% Metric 2B: Lead/Lag via Centroid (Median) Difference
            results.centroidMedianDiff(c, h) = eyeT50 - stimT50;
            
            %% Metric 2C: Temporal Skew from Centroid (Median)
            results.eye.skewFromCentroidMedian(c, h) = 1 - (2 * eyeT50 / halfDuration);
            results.stim.skewFromCentroidMedian(c, h) = 1 - (2 * stimT50 / halfDuration);
            
            %% Metric 2E: Quantile Skew (Bowley Skewness)
            eyeIQR = eyeT75 - eyeT25;
            stimIQR = stimT75 - stimT25;
            
            if eyeIQR > 0
                results.eye.skewQuantile(c, h) = ((eyeT75 - eyeT50) - (eyeT50 - eyeT25)) / eyeIQR;
            else
                results.eye.skewQuantile(c, h) = NaN;
            end
            
            if stimIQR > 0
                results.stim.skewQuantile(c, h) = ((stimT75 - stimT50) - (stimT50 - stimT25)) / stimIQR;
            else
                results.stim.skewQuantile(c, h) = NaN;
            end
            
            %% Metric 3A: Mean-Median Skew
            if eyeIQR > 0
                results.eye.skewMeanMedian(c, h) = (eyeCentroidMean - eyeT50) / eyeIQR;
            else
                results.eye.skewMeanMedian(c, h) = NaN;
            end
            
            if stimIQR > 0
                results.stim.skewMeanMedian(c, h) = (stimCentroidMean - stimT50) / stimIQR;
            else
                results.stim.skewMeanMedian(c, h) = NaN;
            end

            %% Metric 4A: Two Tallest Peaks
            if isequal(h,1), sgn=1; elseif isequal(h,2), sgn=-1; end
            [eyePks, eyePkTimes] = findpeaks(sgn*eyeHalf, tHalf, ...
                'SortStr','descend', 'NPeaks',2, ...
                'MinPeakHeight',0.01, 'MinPeakDistance',0.05);
            if numel(eyePks) >= 1
                results.eye.peak1Val(c, h)  = sgn*eyePks(1);
                results.eye.peak1TimeMs(c, h) = eyePkTimes(1);
            end
            if numel(eyePks) >= 2
                results.eye.peak2Val(c, h)  = sgn*eyePks(2);
                results.eye.peak2TimeMs(c, h) = eyePkTimes(2);
            end
        end
    end
    
    %% Add metadata
    results.metadata.nCycles = nCycles;
    results.metadata.cycleLength = cycleLength;
    results.metadata.halfCycleLength = halfLength;
    results.metadata.units.time = 'ms (assuming 1000 Hz sample rate)';
    results.metadata.units.skewIndices = 'dimensionless';
    results.metadata.columnLabels = {'HalfCycle1', 'HalfCycle2'};
    results.metadata.signConventions.centroidDiff = 'Positive = eye lags stimulus, Negative = eye leads stimulus';
    results.metadata.signConventions.skewFromCentroid = 'Positive = mass early, Negative = mass late';
    results.metadata.signConventions.skewQuantile = 'Positive = late tail, Negative = early tail';
    results.metadata.signConventions.skewMeanMedian = 'Positive = mean later than median, Negative = mean earlier than median';
    results.metadata.description = ["", "CENTROID (MEAN) BASED METRICS:", ...
        "    - These metrics use the temporal center of mass, which is analogous to the mean of a probability distribution. Sensitive to extreme values and tails.", ...
        "    - CENTROID: Eye velocity center-of-mass time (ms). Earlier values indicate mass concentrated toward the start of the half-cycle; later values indicate mass concentrated toward the end.", ...
        "    - LAG via Centroid Difference: Measure of the temporal shift between eye response and stimulus centers-of-mass. Positive values indicate the eye response lags the stimulus; negative values indicate the eye response leads.", ...
        "    - SKEW via Centroid: Measures where the center-of-mass falls relative to the temporal midpoint of the half-cycle. Positive or negative values correspond to early or late mass concentration.", ...
        "", "CENTROID (MEDIAN) AND QUANTILE BASED METRICS:", ...
        "    - These metrics use cumulative area to find quantile times, analogous to percentiles of a probability distribution. Robust to extreme values.", ...
        "    - CENTROID: Eye velocity 50th percentile time (ms). The time at which 50% of the total response area has occurred. Unlike the mean-based centroid, this is robust to outlier peaks or long tails.", ...
        "    - LAG via Centroid Difference: Same as above but using the median-based centroid. More robust to waveform distortions. Positive = lag, Negative = lead.", ...
        "    - SKEW via Centroid: Same as above but using the median-based centroid. Positive = median early, Negative = median late.", ...
        "    - BOWLEY SKEWNESS (Quantile Skew): Measures asymmetry in the spread around the median. Positive = late/right tail, Negative = early/left tail.", ...
        "", "MEAN-MEDIAN COMPARISON METRIC:", ...
        "    - Comparing centroid (mean) to centroid (median) reveals the influence of tails and outliers on the waveform's central tendency.", ...
        "    - MEAN-MEDIAN SKEW: Quantifies how much tails pull the mean away from the median, normalized by the spread of the distribution. Positive = mean is later than median, Negative = mean is earlier than median."];
end

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function [data, nCycles, cycleLength] = orientData(data)
% ORIENTDATA Ensures data is oriented as CYCLELENGTHxNCYCLES
%
% Assumes cycle length is always the larger dimension.

    [nRows, nCols] = size(data);
    
    if nRows < nCols
        % Data is NCYCLESxCYCLELENGTH, transpose it
        data = data';
        cycleLength = nCols;
        nCycles = nRows;
    else
        % Data is already CYCLELENGTHxNCYCLES
        cycleLength = nRows;
        nCycles = nCols;
    end
end

function centroid = computeCentroid(signal, t)
% COMPUTECENTROID Computes temporal center of mass (mean)
%
%   centroid = sum(signal .* t) / sum(signal)
%
% This is analogous to the expected value (mean) of a probability distribution
% where the signal represents the probability density.

    signalSum = sum(signal);
    if signalSum > 0
        centroid = sum(signal .* t) / signalSum;
    else
        centroid = NaN;
    end
end

function [t25, t50, t75] = computeQuantileTimes(signal, t)
% COMPUTEQUANTILETIMES Finds times at which cumulative area reaches 25%, 50%, 75%
%
% Uses the trapezoidal rule for integration and linear interpolation for
% sub-sample precision.
%
% t50 (the median) is analogous to the centroid but uses the median instead
% of the mean, making it robust to outliers and extreme values.

    % Compute cumulative integral using trapezoidal rule
    cumArea = cumtrapz(t, signal);
    
    % Normalize to [0, 1]
    totalArea = cumArea(end);
    
    if totalArea <= 0
        t25 = NaN;
        t50 = NaN;
        t75 = NaN;
        return;
    end
    
    cumAreaNorm = cumArea / totalArea;
    
    % Find quantile times via interpolation
    t25 = interp1(cumAreaNorm, t, 0.25, 'linear');
    t50 = interp1(cumAreaNorm, t, 0.50, 'linear');
    t75 = interp1(cumAreaNorm, t, 0.75, 'linear');
end