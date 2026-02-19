function results = cycleVarianceDecomposition(cycles, varargin)
% CYCLEVARIANCEDECOMPOSITION  Variance decomposition of sinusoidal cycle data.
%
%   results = CYCLEVARIANCEDECOMPOSITION(cycles)
%   results = CYCLEVARIANCEDECOMPOSITION(cycles, 'Name', Value, ...)
%
%   Decomposes the pointwise variance of a matrix of individual cycles into
%   sinusoidal (fundamental) and residual (non-sinusoidal) components. This
%   is useful for characterizing the error structure in cycle-averaged
%   vestibulo-ocular reflex (VOR) data or any periodic signal.
%
%   The function is NaN-aware throughout: NaN values (from blink removal,
%   saccade excision, tracker dropout, etc.) are handled gracefully at
%   every stage. Per-cycle sinusoidal fits use only valid (non-NaN) samples,
%   pointwise statistics are computed with NaN-omission, and the effective
%   sample count is tracked at every time point and per cycle.
%
%   The function automatically detects the orientation of the input matrix
%   by assuming the longer dimension corresponds to the number of cycles
%   (since cycle length is typically <= 2000 samples for physiological
%   recordings, while the number of cycles is often >= 30). If the matrix
%   is square or the assumption is ambiguous, use the 'Orientation' option.
%
%   INPUTS
%   ------
%   cycles : [nCycles x cycleLength] or [cycleLength x nCycles] matrix
%       Matrix of individual cycles. Each cycle is one period of the
%       expected sinusoidal response. NaN values are permitted.
%
%   NAME-VALUE PARAMETERS
%   ---------------------
%   'Fs'            - Sampling rate in Hz (default: 1000)
%   'StimFreq'      - Stimulus frequency in Hz (default: 1)
%   'Orientation'   - 'auto' (default), 'rows', or 'cols'
%                     'rows' = each row is a cycle
%                     'cols' = each column is a cycle
%                     'auto' = longer dimension is assumed to be nCycles
%   'Alpha'         - Significance level for confidence intervals
%                     (default: 0.05 for 95% CI)
%   'MinValidFrac'  - Minimum fraction of valid (non-NaN) samples required
%                     in a cycle for its sinusoidal fit to be attempted
%                     (default: 0.5). Cycles below this threshold are
%                     excluded from coefficient-level analyses.
%   'MinValidCycles'- Minimum number of valid cycles at a time point for
%                     pointwise statistics to be computed (default: 3).
%                     Time points with fewer valid cycles get NaN.
%
%   OUTPUT
%   ------
%   results : struct with fields organized into logical groups:
%
%     --- Input metadata ---
%     .nCycles          : Total number of cycles in input
%     .cycleLength      : Number of samples per cycle
%     .Fs               : Sampling rate (Hz)
%     .stimFreq         : Stimulus frequency (Hz)
%     .alpha            : Significance level used for CIs
%     .orientation      : Detected or specified orientation
%     .timeAxis         : [1 x cycleLength] time vector (seconds)
%
%     --- NaN diagnostics ---
%     .nValidCyclesAt   : [1 x cycleLength] valid cycle count per timepoint
%     .nValidSamplesIn  : [nCycles x 1] valid sample count per cycle
%     .fracValidIn      : [nCycles x 1] fraction of valid samples per cycle
%     .cycleIncluded    : [nCycles x 1] logical, true if cycle met the
%                         MinValidFrac threshold for sinusoidal fitting
%     .nCyclesIncluded  : Number of cycles included in coefficient analyses
%     .nCyclesExcluded  : Number of cycles excluded
%
%     --- Cycle-mean waveform ---
%     .cycleMean        : [1 x cycleLength] mean waveform (NaN-omitting)
%     .cycleSEM         : [1 x cycleLength] pointwise SEM (total)
%     .cycleMeanCI      : [2 x cycleLength] lower and upper CI bounds
%
%     --- Sinusoidal fit of the cycle-mean ---
%     .fitMean          : [1 x cycleLength] sinusoidal fit of cycle-mean
%     .fitMeanCoeffs    : struct with fields A, B, C (sin, cos, offset)
%     .fitMeanAmplitude : Scalar amplitude of cycle-mean fit
%     .fitMeanPhase_deg : Scalar phase of cycle-mean fit (degrees)
%
%     --- Per-cycle sinusoidal fits (all nCycles; excluded = NaN) ---
%     .fitCoeffsAll     : [nCycles x 3] matrix of [A, B, C] per cycle
%     .fitAmplitudeAll  : [nCycles x 1] amplitude per cycle
%     .fitPhaseAll_deg  : [nCycles x 1] phase per cycle (degrees)
%     .fitWaveformsAll  : [nCycles x cycleLength] fitted sinusoid per cycle
%     .residualsAll     : [nCycles x cycleLength] residuals per cycle
%     .fitRsquaredAll   : [nCycles x 1] R-squared of each cycle's fit
%
%     --- Coefficient statistics (included cycles only) ---
%     .coeffMean        : [1 x 3] mean of [A, B, C] across included cycles
%     .coeffSEM         : [1 x 3] SEM of [A, B, C] across included cycles
%     .coeffCov         : [3 x 3] covariance matrix of [A, B, C]
%
%     --- Amplitude and phase uncertainty (delta method) ---
%     .amplitudeMean    : Scalar amplitude from mean coefficients
%     .amplitudeSEM     : Scalar SEM of amplitude (delta method)
%     .amplitudeCI      : [1 x 2] CI for amplitude
%     .phaseMean_deg    : Scalar phase from mean coefficients (degrees)
%     .phaseSEM_deg     : Scalar SEM of phase (delta method, degrees)
%     .phaseCI_deg      : [1 x 2] CI for phase (degrees)
%
%     --- Variance decomposition (pointwise, 1 x cycleLength arrays) ---
%     .varTotal         : Total pointwise variance across cycles
%     .varSinusoidal    : Pointwise variance of sinusoidal components
%     .varResidual      : Pointwise variance of residual components
%     .varCrossterm     : Pointwise cross-term (2*Cov[s,e])
%     .semSinusoidal    : Pointwise SEM of sinusoidal component
%     .semResidual      : Pointwise SEM of residual component
%     .fractSinusoidal  : Pointwise fraction of variance (sinusoidal)
%     .fractResidual    : Pointwise fraction of variance (residual)
%
%     --- Global variance summary (scalars, averaged over time) ---
%     .varTotalMean         : Mean total variance across cycle
%     .varSinusoidalMean    : Mean sinusoidal variance across cycle
%     .varResidualMean      : Mean residual variance across cycle
%     .etaSquared           : Overall fraction of variance from sinusoidal
%     .residualFraction     : Overall fraction of variance from residual
%
%     --- Data quality diagnostics ---
%     .snr_dB               : Pointwise SNR in dB (sinusoidal/residual var)
%     .snrMean_dB           : Mean SNR across cycle (dB)
%     .residualRMS          : [nCycles x 1] RMS of residuals per cycle
%     .residualRMSMean      : Mean residual RMS across included cycles
%     .residualRMSOutlierTh : Outlier threshold (median + 3*MAD)
%     .outlierCycleIdx      : Indices of suspected outlier cycles
%
%   EXAMPLE
%   -------
%     % Simulate 60 noisy VOR cycles at 1 Hz, Fs = 1000, with NaN gaps
%     Fs = 1000; f = 1; nCycles = 60;
%     cycleLen = Fs / f;
%     t = (0:cycleLen-1) / Fs;
%     cycles = zeros(nCycles, cycleLen);
%     for k = 1:nCycles
%         amp = 50 + 5*randn;
%         phi = 0.1*randn;
%         cycles(k,:) = amp*sin(2*pi*f*t + phi) + 5*randn(1,cycleLen);
%     end
%     % Simulate blink artifacts (NaN gaps) in ~30% of cycles
%     for k = 1:nCycles
%         if rand < 0.3
%             blinkStart = randi([100, 800]);
%             blinkLen = randi([50, 150]);
%             blinkEnd = min(blinkStart + blinkLen - 1, cycleLen);
%             cycles(k, blinkStart:blinkEnd) = NaN;
%         end
%     end
%     results = cycleVarianceDecomposition(cycles, 'Fs', 1000, 'StimFreq', 1);
%
%     % Plot cycle-mean with CI band and valid-count overlay
%     figure;
%     subplot(2,1,1);
%     fill([results.timeAxis fliplr(results.timeAxis)], ...
%          [results.cycleMeanCI(1,:) fliplr(results.cycleMeanCI(2,:))], ...
%          [0.8 0.8 1], 'EdgeColor', 'none'); hold on;
%     plot(results.timeAxis, results.cycleMean, 'b', 'LineWidth', 1.5);
%     plot(results.timeAxis, results.fitMean, 'r--', 'LineWidth', 1.2);
%     xlabel('Time (s)'); ylabel('Eye velocity (deg/s)');
%     legend('95% CI', 'Cycle mean', 'Sinusoidal fit');
%     title(sprintf('Cycle Mean (%.0f included / %.0f total cycles)', ...
%         results.nCyclesIncluded, results.nCycles));
%     subplot(2,1,2);
%     plot(results.timeAxis, results.nValidCyclesAt, 'k', 'LineWidth',1);
%     xlabel('Time (s)'); ylabel('Valid cycles');
%     title('Valid cycle count per time point');
%
%   See also MEAN, STD, VAR, COV, ATAN2
%
%   Author:  Generated with assistance from Claude (Anthropic)
%   Version: 2.0 (NaN-aware)

    %% Parse inputs
    p = inputParser;
    addRequired(p, 'cycles', @(x) isnumeric(x) && ismatrix(x));
    addParameter(p, 'Fs', 1000, @(x) isscalar(x) && x > 0);
    addParameter(p, 'StimFreq', 1, @(x) isscalar(x) && x > 0);
    addParameter(p, 'Orientation', 'auto', ...
        @(x) ismember(lower(x), {'auto','rows','cols'}));
    addParameter(p, 'Alpha', 0.05, @(x) isscalar(x) && x > 0 && x < 1);
    addParameter(p, 'MinValidFrac', 0.5, ...
        @(x) isscalar(x) && x > 0 && x <= 1);
    addParameter(p, 'MinValidCycles', 3, ...
        @(x) isscalar(x) && x >= 2 && round(x) == x);
    parse(p, cycles, varargin{:});

    Fs              = p.Results.Fs;
    stimFreq        = p.Results.StimFreq;
    orient          = lower(p.Results.Orientation);
    alpha           = p.Results.Alpha;
    minValidFrac    = p.Results.MinValidFrac;
    minValidCycles  = p.Results.MinValidCycles;

    %% Orient matrix to [nCycles x cycleLength]
    [nRows, nCols] = size(cycles);

    switch orient
        case 'auto'
            expectedCycleLen = round(Fs / stimFreq);
            if nCols == expectedCycleLen && nRows ~= expectedCycleLen
                orientUsed = 'rows';
            elseif nRows == expectedCycleLen && nCols ~= expectedCycleLen
                cycles = cycles';
                orientUsed = 'cols';
            elseif nRows == nCols
                warning('cycle_decomposition:squareMatrix', ...
                    ['Square matrix detected (%d x %d). Assuming rows = cycles. ' ...
                     'Use ''Orientation'' parameter if this is incorrect.'], ...
                    nRows, nCols);
                orientUsed = 'rows (assumed, square matrix)';
            else
                if nRows > nCols
                    orientUsed = 'cols (auto, longer dim)';
                else
                    cycles = cycles';
                    orientUsed = 'rows (auto, longer dim)';
                end
            end
        case 'rows'
            orientUsed = 'rows';
        case 'cols'
            cycles = cycles';
            orientUsed = 'cols';
    end

    [nCycles, cycleLength] = size(cycles);

    % Validate cycle length against expected
    expectedCycleLen = round(Fs / stimFreq);
    if cycleLength ~= expectedCycleLen
        warning('cycle_decomposition:cycleLengthMismatch', ...
            ['Cycle length (%d) does not match expected length (%d) from ' ...
             'Fs=%.1f and StimFreq=%.1f. Proceeding anyway.'], ...
            cycleLength, expectedCycleLen, Fs, stimFreq);
    end

    %% Build time vector and design matrix
    t = (0:cycleLength - 1)' / Fs;           % [cycleLength x 1]
    omega = 2 * pi * stimFreq;
    X = [sin(omega * t), cos(omega * t), ones(cycleLength, 1)];  % [cycleLength x 3]

    %% NaN diagnostics
    validMask       = ~isnan(cycles);                             % [nCycles x cycleLength]
    nValidCyclesAt  = sum(validMask, 1);                          % [1 x cycleLength]
    nValidSamplesIn = sum(validMask, 2);                          % [nCycles x 1]
    fracValidIn     = nValidSamplesIn / cycleLength;              % [nCycles x 1]
    cycleIncluded   = fracValidIn >= minValidFrac;                % [nCycles x 1] logical
    nCyclesIncluded = sum(cycleIncluded);
    nCyclesExcluded = nCycles - nCyclesIncluded;

    totalNaN = sum(~validMask(:));
    totalElem = numel(cycles);
    if totalNaN > 0
        fprintf(['  cycle_decomposition: %.1f%% of all samples are NaN ' ...
                 '(%d/%d).\n'], 100*totalNaN/totalElem, totalNaN, totalElem);
        fprintf('  %d/%d cycles meet MinValidFrac threshold (%.0f%%).\n', ...
            nCyclesIncluded, nCycles, minValidFrac*100);
    end

    if nCyclesIncluded < minValidCycles
        warning('cycle_decomposition:tooFewCycles', ...
            ['Only %d cycles meet the MinValidFrac threshold. ' ...
             'Coefficient-level analyses require at least %d. ' ...
             'Consider lowering MinValidFrac or checking data quality.'], ...
            nCyclesIncluded, minValidCycles);
    end

    %% Cycle-mean and pointwise SEM (NaN-omitting)
    cycleMean = mean(cycles, 1, 'omitnan');                       % [1 x cycleLength]
    cycleSD   = std(cycles, 0, 1, 'omitnan');                     % [1 x cycleLength]
    cycleSEM  = cycleSD ./ sqrt(max(nValidCyclesAt, 1));           % [1 x cycleLength]

    % Apply MinValidCycles threshold: mask time points with too few cycles
    tooFewMask = nValidCyclesAt < minValidCycles;
    cycleMean(tooFewMask) = NaN;
    cycleSEM(tooFewMask)  = NaN;

    % CI uses t-distribution with (nValid - 1) degrees of freedom per point
    tCritVec = zeros(1, cycleLength);
    for j = 1:cycleLength
        if nValidCyclesAt(j) >= minValidCycles
            tCritVec(j) = tinv(1 - alpha/2, nValidCyclesAt(j) - 1);
        else
            tCritVec(j) = NaN;
        end
    end
    cycleMeanCI = [cycleMean - tCritVec .* cycleSEM; ...
                   cycleMean + tCritVec .* cycleSEM];              % [2 x cycleLength]

    %% Sinusoidal fit of the cycle-mean (using valid time points only)
    validMeanMask = ~isnan(cycleMean);
    if sum(validMeanMask) >= 3  % need at least 3 points for 3 parameters
        Xm = X(validMeanMask', :);
        ym = cycleMean(validMeanMask)';
        betaMean = Xm \ ym;
        fitMean  = (X * betaMean)';                                % [1 x cycleLength]
    else
        betaMean = [NaN; NaN; NaN];
        fitMean  = NaN(1, cycleLength);
        warning('cycle_decomposition:noValidMean', ...
            'Too few valid time points in cycle-mean for sinusoidal fit.');
    end

    fitMeanCoeffs.A = betaMean(1);
    fitMeanCoeffs.B = betaMean(2);
    fitMeanCoeffs.C = betaMean(3);
    fitMeanAmplitude = sqrt(betaMean(1)^2 + betaMean(2)^2);
    fitMeanPhase_deg = atan2d(betaMean(2), betaMean(1));

    %% Per-cycle sinusoidal fits (NaN-aware)
    fitCoeffsAll    = NaN(nCycles, 3);
    fitWaveformsAll = NaN(nCycles, cycleLength);
    residualsAll    = NaN(nCycles, cycleLength);
    fitRsquaredAll  = NaN(nCycles, 1);

    for k = 1:nCycles
        validK = validMask(k, :)';                                % [cycleLength x 1] logical
        nValidK = sum(validK);

        if ~cycleIncluded(k) || nValidK < 3
            % Skip this cycle — output remains NaN
            continue;
        end

        Xk = X(validK, :);
        yk = cycles(k, validK)';

        beta_k = Xk \ yk;
        fitCoeffsAll(k, :) = beta_k';

        % Reconstruct full-length fitted waveform (defined at all points)
        fitWaveformsAll(k, :) = (X * beta_k)';

        % Residuals only at valid sample positions
        residualsAll(k, validK) = cycles(k, validK) - fitWaveformsAll(k, validK);

        % R-squared for this cycle (computed on valid samples only)
        ssRes = sum((yk - Xk * beta_k).^2);
        ssTot = sum((yk - mean(yk)).^2);
        if ssTot > 0
            fitRsquaredAll(k) = 1 - ssRes / ssTot;
        else
            fitRsquaredAll(k) = NaN;
        end
    end

    fitAmplitudeAll = sqrt(fitCoeffsAll(:,1).^2 + fitCoeffsAll(:,2).^2);
    fitPhaseAll_deg = atan2d(fitCoeffsAll(:,2), fitCoeffsAll(:,1));

    %% Coefficient statistics across included cycles
    inclCoeffs = fitCoeffsAll(cycleIncluded, :);                   % [nIncl x 3]

    if nCyclesIncluded >= minValidCycles
        coeffMean = mean(inclCoeffs, 1, 'omitnan');
        coeffSEM  = std(inclCoeffs, 0, 1, 'omitnan') / sqrt(nCyclesIncluded);
        coeffCov  = nancov(inclCoeffs);
    else
        coeffMean = NaN(1, 3);
        coeffSEM  = NaN(1, 3);
        coeffCov  = NaN(3, 3);
    end

    %% Amplitude and phase uncertainty via delta method
    A_bar = coeffMean(1);
    B_bar = coeffMean(2);
    ampMean = sqrt(A_bar^2 + B_bar^2);

    if nCyclesIncluded >= minValidCycles && ampMean > 0
        covAB = coeffCov(1:2, 1:2);

        % Gradient of amplitude w.r.t. [A, B]
        gradAmp = [A_bar; B_bar] / ampMean;
        varAmp  = gradAmp' * (covAB / nCyclesIncluded) * gradAmp;
        semAmp  = sqrt(max(varAmp, 0));

        % Gradient of phase = atan2(B, A) w.r.t. [A, B]
        gradPhi    = [-B_bar; A_bar] / ampMean^2;
        varPhi     = gradPhi' * (covAB / nCyclesIncluded) * gradPhi;
        semPhi_rad = sqrt(max(varPhi, 0));
        semPhi_deg = rad2deg(semPhi_rad);

        tCritCoeff  = tinv(1 - alpha/2, nCyclesIncluded - 1);
    else
        semAmp      = NaN;
        semPhi_deg  = NaN;
        tCritCoeff  = NaN;
    end

    phaseMean_deg = atan2d(B_bar, A_bar);
    amplitudeCI   = ampMean       + tCritCoeff * semAmp     * [-1, 1];
    phaseCI_deg   = phaseMean_deg + tCritCoeff * semPhi_deg * [-1, 1];

    %% Variance decomposition (pointwise, using included cycles only)
    inclFitWaveforms = fitWaveformsAll(cycleIncluded, :);
    inclResiduals    = residualsAll(cycleIncluded, :);

    % For pointwise variance, we need to handle NaN in the residuals
    % (residuals are NaN where original data was NaN)
    varTotal      = var(cycles(cycleIncluded, :), 0, 1, 'omitnan');
    varSinusoidal = var(inclFitWaveforms, 0, 1, 'omitnan');
    varResidual   = var(inclResiduals, 0, 1, 'omitnan');

    varCrossterm  = varTotal - varSinusoidal - varResidual;

    % Pointwise effective N for SEM calculation (valid residuals at each t)
    nValidResidAt = sum(~isnan(inclResiduals), 1);
    nValidFitAt   = sum(~isnan(inclFitWaveforms), 1);

    semSinusoidal = sqrt(varSinusoidal) ./ sqrt(max(nValidFitAt, 1));
    semResidual   = sqrt(varResidual)   ./ sqrt(max(nValidResidAt, 1));

    % Fraction of variance (guard against division by zero)
    safeVarTotal    = max(varTotal, eps);
    fractSinusoidal = varSinusoidal ./ safeVarTotal;
    fractResidual   = varResidual   ./ safeVarTotal;

    % Apply MinValidCycles mask to decomposition outputs
    decompTooFew = nValidResidAt < minValidCycles;
    varTotal(decompTooFew)         = NaN;
    varSinusoidal(decompTooFew)    = NaN;
    varResidual(decompTooFew)      = NaN;
    varCrossterm(decompTooFew)     = NaN;
    semSinusoidal(decompTooFew)    = NaN;
    semResidual(decompTooFew)      = NaN;
    fractSinusoidal(decompTooFew)  = NaN;
    fractResidual(decompTooFew)    = NaN;

    %% Global variance summary (scalars, averaged over valid time points)
    varTotalMean      = mean(varTotal, 'omitnan');
    varSinusoidalMean = mean(varSinusoidal, 'omitnan');
    varResidualMean   = mean(varResidual, 'omitnan');
    etaSquared        = varSinusoidalMean / max(varTotalMean, eps);
    residualFraction  = varResidualMean   / max(varTotalMean, eps);

    %% Data quality diagnostics
    % Pointwise SNR
    snr_dB     = 10 * log10(varSinusoidal ./ max(varResidual, eps));
    snrMean_dB = 10 * log10(varSinusoidalMean / max(varResidualMean, eps));

    % Per-cycle residual RMS for outlier detection (NaN-aware)
    residualRMS = sqrt(mean(residualsAll.^2, 2, 'omitnan'));      % [nCycles x 1]

    inclResidRMS = residualRMS(cycleIncluded);
    inclResidRMS = inclResidRMS(~isnan(inclResidRMS));

    if numel(inclResidRMS) >= 3
        residualRMSMean = mean(inclResidRMS);
        medRMS    = median(inclResidRMS);
        madRMS    = median(abs(inclResidRMS - medRMS));
        outlierTh = medRMS + 3 * 1.4826 * madRMS;
        outlierIdx = find(residualRMS > outlierTh & cycleIncluded);
    else
        residualRMSMean = NaN;
        outlierTh       = NaN;
        outlierIdx      = [];
    end

    %% Pack results struct
    results = struct();

    % --- Input metadata ---
    results.nCycles          = nCycles;
    results.cycleLength      = cycleLength;
    results.Fs               = Fs;
    results.stimFreq         = stimFreq;
    results.alpha            = alpha;
    results.orientation      = orientUsed;
    results.timeAxis         = t';

    % --- NaN diagnostics ---
    results.nValidCyclesAt   = nValidCyclesAt;
    results.nValidSamplesIn  = nValidSamplesIn;
    results.fracValidIn      = fracValidIn;
    results.cycleIncluded    = cycleIncluded;
    results.nCyclesIncluded  = nCyclesIncluded;
    results.nCyclesExcluded  = nCyclesExcluded;

    % --- Cycle-mean waveform ---
    results.cycleMean        = cycleMean;
    results.cycleSEM         = cycleSEM;
    results.cycleMeanCI      = cycleMeanCI;

    % --- Sinusoidal fit of the cycle-mean ---
    results.fitMean          = fitMean;
    results.fitMeanCoeffs    = fitMeanCoeffs;
    results.fitMeanAmplitude = fitMeanAmplitude;
    results.fitMeanPhase_deg = fitMeanPhase_deg;

    % --- Per-cycle sinusoidal fits ---
    results.fitCoeffsAll     = fitCoeffsAll;
    results.fitAmplitudeAll  = fitAmplitudeAll;
    results.fitPhaseAll_deg  = fitPhaseAll_deg;
    results.fitWaveformsAll  = fitWaveformsAll;
    results.residualsAll     = residualsAll;
    results.fitRsquaredAll   = fitRsquaredAll;

    % --- Coefficient statistics ---
    results.coeffMean        = coeffMean;
    results.coeffSEM         = coeffSEM;
    results.coeffCov         = coeffCov;

    % --- Amplitude and phase uncertainty ---
    results.amplitudeMean    = ampMean;
    results.amplitudeSEM     = semAmp;
    results.amplitudeCI      = amplitudeCI;
    results.phaseMean_deg    = phaseMean_deg;
    results.phaseSEM_deg     = semPhi_deg;
    results.phaseCI_deg      = phaseCI_deg;

    % --- Variance decomposition (pointwise) ---
    results.varTotal         = varTotal;
    results.varSinusoidal    = varSinusoidal;
    results.varResidual      = varResidual;
    results.varCrossterm     = varCrossterm;
    results.semSinusoidal    = semSinusoidal;
    results.semResidual      = semResidual;
    results.fractSinusoidal  = fractSinusoidal;
    results.fractResidual    = fractResidual;

    % --- Global variance summary ---
    results.varTotalMean         = varTotalMean;
    results.varSinusoidalMean    = varSinusoidalMean;
    results.varResidualMean      = varResidualMean;
    results.etaSquared           = etaSquared;
    results.residualFraction     = residualFraction;

    % --- Data quality diagnostics ---
    results.snr_dB               = snr_dB;
    results.snrMean_dB           = snrMean_dB;
    results.residualRMS          = residualRMS;
    results.residualRMSMean      = residualRMSMean;
    results.residualRMSOutlierTh = outlierTh;
    results.outlierCycleIdx      = outlierIdx;

end

%% Local helper: NaN-aware covariance
function C = nancov(X)
% NANCOV  Covariance matrix with listwise NaN deletion.
%   Removes any row of X containing a NaN, then computes cov().
    validRows = all(~isnan(X), 2);
    if sum(validRows) < 2
        C = NaN(size(X, 2));
    else
        C = cov(X(validRows, :));
    end
end