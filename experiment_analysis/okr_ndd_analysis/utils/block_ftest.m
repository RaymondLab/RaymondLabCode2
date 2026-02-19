function results = block_ftest(cycles, blockLabel, varargin)
% BLOCK_FTEST  Pointwise omnibus F-test for block effects across cycles.
%
%   results = BLOCK_FTEST(cycles, blockLabel)
%   results = BLOCK_FTEST(cycles, blockLabel, 'Name', Value, ...)
%
%   Performs a one-way ANOVA at each time point across N blocks to test
%   whether blocks are exchangeable replications of the same response.
%   If the omnibus test is significant, post-hoc pairwise t-tests identify
%   which blocks differ. NaN values are handled gracefully throughout.
%
%   The default multiple comparisons correction is cluster-based
%   permutation testing, which respects the temporal correlation structure
%   of the data and is sensitive to spatially extended effects (e.g.,
%   sub-regions where block means diverge). Pointwise corrections
%   (Bonferroni, FDR) are also available but are more conservative for
%   temporally correlated data.
%
%   INPUTS
%   ------
%   cycles     : [nCycles x cycleLength] matrix of individual cycles.
%                NaN values are permitted.
%   blockLabel : [nCycles x 1] vector of block identifiers (numeric or
%                categorical). Each unique value defines a block.
%
%   NAME-VALUE PARAMETERS
%   ---------------------
%   'Alpha'          - Significance level (default: 0.05)
%   'Correction'     - Multiple comparisons correction for pointwise tests:
%                      'cluster'    : Cluster-based permutation (default)
%                      'none'       : No correction
%                      'bonferroni' : Bonferroni correction
%                      'fdr'        : Benjamini-Hochberg FDR
%                      (default: 'cluster')
%   'MinValidPerBlock' - Minimum valid (non-NaN) cycles required per block
%                        at a time point for the test to be computed
%                        (default: 3). Time points where any block has
%                        fewer valid cycles receive NaN.
%   'Fs'             - Sampling rate in Hz, used to construct time axis
%                      (default: 1000)
%   'StimFreq'       - Stimulus frequency in Hz (default: 1)
%   'PostHoc'        - Post-hoc method when omnibus is significant:
%                      'tukey'  : Tukey's HSD (balanced designs only)
%                      'holm'   : Holm-Bonferroni corrected t-tests
%                      (default: 'holm')
%   'nPermutations'  - Number of permutations for cluster correction
%                      (default: 5000). Ignored for other corrections.
%   'ClusterAlpha'   - Cluster-forming threshold: the uncorrected alpha
%                      level used to threshold the F-statistics before
%                      forming clusters (default: 0.05). Lower values
%                      form fewer, more conservative clusters. Ignored
%                      for non-cluster corrections.
%
%   OUTPUT
%   ------
%   results : struct with fields:
%
%     --- Metadata ---
%     .nBlocks            : Number of unique blocks
%     .blockIDs           : [nBlocks x 1] unique block identifiers
%     .nCyclesPerBlock    : [nBlocks x 1] total cycles per block
%     .cycleLength        : Number of samples per cycle
%     .alpha              : Significance level used
%     .correction         : Multiple comparisons method used
%     .postHocMethod      : Post-hoc method used
%     .timeAxis           : [1 x cycleLength] time vector (seconds)
%
%     --- Block-level descriptives ---
%     .blockMeans         : [nBlocks x cycleLength] cycle-mean per block
%     .blockSEMs          : [nBlocks x cycleLength] SEM per block
%     .blockNValid        : [nBlocks x cycleLength] valid cycles per block
%                           at each time point
%     .grandMean          : [1 x cycleLength] overall cycle-mean
%
%     --- Omnibus F-test (pointwise) ---
%     .fStats             : [1 x cycleLength] F-statistic at each point
%     .pValues            : [1 x cycleLength] raw (uncorrected) p-values
%     .pCorrected         : [1 x cycleLength] corrected p-values
%                           (NaN for cluster correction — see sigClusters)
%     .significant        : [1 x cycleLength] logical, true if significant
%                           after correction. For cluster correction, time
%                           points belonging to any significant cluster.
%     .dfBetween          : Numerator degrees of freedom (nBlocks - 1)
%     .dfWithinAt         : [1 x cycleLength] denominator df per timepoint
%     .etaSquaredPartial  : [1 x cycleLength] partial eta-squared effect
%                           size (SS_between / SS_total) at each point
%
%     --- Summary statistics ---
%     .nSigPoints         : Number of significant time points
%     .fracSigPoints      : Fraction of cycle that is significant
%     .maxF               : Maximum F-statistic across cycle
%     .maxFTime           : Time (s) at which maximum F occurs
%     .sigClusters        : Struct array of significant clusters with:
%                           .startIdx, .endIdx    - sample indices
%                           .startTime, .endTime  - in seconds
%                           .duration             - cluster duration (s)
%                           .meanF, .peakF        - F-statistic summaries
%                           .mass                 - sum of F in cluster
%                           .pValue               - cluster-level p-value
%                                                   (only for 'cluster')
%     .allClusters        : Struct array of ALL supra-threshold clusters,
%                           including non-significant (cluster only)
%     .recommendation     : Text string summarizing the result
%
%     --- Cluster permutation details (only for 'cluster' correction) ---
%     .clusterAlpha       : Cluster-forming alpha used
%     .clusterFThreshold  : F-statistic threshold for cluster formation
%     .nPermutations      : Number of permutations performed
%     .nullMaxClusterMass : [nPermutations x 1] null distribution of
%                           maximum cluster masses
%
%     --- Post-hoc pairwise comparisons ---
%     .pairwise           : Struct with fields (only populated at time
%                           points where omnibus is significant):
%       .pairs            : [nPairs x 2] matrix of block ID pairs
%       .tStats           : [nPairs x cycleLength] t-statistics
%       .pValues          : [nPairs x cycleLength] raw pairwise p-values
%       .pCorrected       : [nPairs x cycleLength] corrected within each
%                           time point using the specified post-hoc method
%       .significant      : [nPairs x cycleLength] logical
%       .meanDiffs        : [nPairs x cycleLength] difference in block
%                           means (block_i - block_j) at each time point
%
%   EXAMPLE
%   -------
%     % Simulate 3 blocks of 60 cycles with sub-regional divergence
%     Fs = 1000; f = 1;
%     cycleLen = Fs / f;
%     t = (0:cycleLen-1) / Fs;
%     nCycPerBlock = 60;
%     cycles = []; blockLabel = [];
%     for b = 1:3
%         for k = 1:nCycPerBlock
%             amp = 50 + 5*randn;
%             cyc = amp*sin(2*pi*f*t) + 5*randn(1,cycleLen);
%             % Add block-specific sub-oscillatory bump in block 3
%             if b == 3
%                 bump = 15 * exp(-0.5*((t - 0.3)/0.05).^2);
%                 cyc = cyc + bump;
%             end
%             cycles = [cycles; cyc];
%             blockLabel = [blockLabel; b];
%         end
%     end
%     results = block_ftest(cycles, blockLabel, 'Fs', 1000, 'StimFreq', 1);
%
%     % Quick visual check
%     figure;
%     subplot(3,1,1); hold on;
%     colors = lines(results.nBlocks);
%     for b = 1:results.nBlocks
%         plot(results.timeAxis, results.blockMeans(b,:), ...
%              'Color', colors(b,:), 'LineWidth', 1.5);
%     end
%     legend(arrayfun(@(x) sprintf('Block %d', x), ...
%         results.blockIDs, 'UniformOutput', false));
%     ylabel('Eye velocity (deg/s)');
%     title('Block cycle-means');
%
%     subplot(3,1,2);
%     plot(results.timeAxis, results.fStats, 'k', 'LineWidth', 1);
%     hold on;
%     sigMask = results.significant;
%     if any(sigMask)
%         plot(results.timeAxis(sigMask), results.fStats(sigMask), ...
%              'r.', 'MarkerSize', 8);
%     end
%     yline(results.clusterFThreshold, 'r--', 'Cluster threshold');
%     xlabel('Time (s)'); ylabel('F-statistic');
%     title(sprintf('Omnibus F-test (%d sig. clusters)', ...
%         length(results.sigClusters)));
%
%     subplot(3,1,3);
%     if ~isempty(results.nullMaxClusterMass)
%         histogram(results.nullMaxClusterMass, 50, ...
%             'FaceColor', [0.7 0.7 0.7]); hold on;
%         for c = 1:length(results.sigClusters)
%             xline(results.sigClusters(c).mass, 'r-', ...
%                 sprintf('Cluster %d (p=%.4f)', c, ...
%                 results.sigClusters(c).pValue), 'LineWidth', 2);
%         end
%         xlabel('Cluster mass'); ylabel('Count');
%         title('Permutation null distribution');
%     end
%
%     disp(results.recommendation);
%
%   See also ANOVA1, TTEST2, MULTCOMPARE
%
%   Author:  Generated with assistance from Claude (Anthropic)
%   Version: 2.0 (cluster-based permutation default)

    %% Parse inputs
    p = inputParser;
    addRequired(p, 'cycles', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'blockLabel', @(x) (isnumeric(x) || iscategorical(x)) ...
        && isvector(x));
    addParameter(p, 'Alpha', 0.05, @(x) isscalar(x) && x > 0 && x < 1);
    addParameter(p, 'Correction', 'cluster', ...
        @(x) ismember(lower(x), {'none','bonferroni','fdr','cluster'}));
    addParameter(p, 'MinValidPerBlock', 3, ...
        @(x) isscalar(x) && x >= 1 && round(x) == x);
    addParameter(p, 'Fs', 1000, @(x) isscalar(x) && x > 0);
    addParameter(p, 'StimFreq', 1, @(x) isscalar(x) && x > 0);
    addParameter(p, 'PostHoc', 'holm', ...
        @(x) ismember(lower(x), {'tukey','holm'}));
    addParameter(p, 'nPermutations', 10000, ...
        @(x) isscalar(x) && x >= 100 && round(x) == x);
    addParameter(p, 'ClusterAlpha', 0.05, ...
        @(x) isscalar(x) && x > 0 && x < 1);
    parse(p, cycles, blockLabel, varargin{:});

    alpha           = p.Results.Alpha;
    correction      = lower(p.Results.Correction);
    minValidPerBlk  = p.Results.MinValidPerBlock;
    Fs              = p.Results.Fs;
    stimFreq        = p.Results.StimFreq;
    postHocMethod   = lower(p.Results.PostHoc);
    nPermutations   = p.Results.nPermutations;
    clusterAlpha    = p.Results.ClusterAlpha;

    %% Validate inputs
    blockLabel = blockLabel(:);  % force column
    if iscategorical(blockLabel)
        blockLabel = double(blockLabel);
    end

    [nCycles, cycleLength] = size(cycles);
    if length(blockLabel) ~= nCycles
        error('block_ftest:dimensionMismatch', ...
            'blockLabel length (%d) must match number of rows in cycles (%d).', ...
            length(blockLabel), nCycles);
    end

    blockIDs = unique(blockLabel);
    nBlocks  = length(blockIDs);

    if nBlocks < 2
        error('block_ftest:insufficientBlocks', ...
            'Need at least 2 blocks for F-test, found %d.', nBlocks);
    end

    nCyclesPerBlock = zeros(nBlocks, 1);
    for b = 1:nBlocks
        nCyclesPerBlock(b) = sum(blockLabel == blockIDs(b));
    end

    timeAxis = (0:cycleLength - 1) / Fs;

    %% Block-level descriptives
    blockMeans  = NaN(nBlocks, cycleLength);
    blockSEMs   = NaN(nBlocks, cycleLength);
    blockNValid = zeros(nBlocks, cycleLength);

    for b = 1:nBlocks
        bCycles = cycles(blockLabel == blockIDs(b), :);
        validMask = ~isnan(bCycles);
        blockNValid(b, :) = sum(validMask, 1);
        blockMeans(b, :)  = mean(bCycles, 1, 'omitnan');
        blockSDs          = std(bCycles, 0, 1, 'omitnan');
        blockSEMs(b, :)   = blockSDs ./ sqrt(max(blockNValid(b,:), 1));
    end

    grandMean = mean(cycles, 1, 'omitnan');

    %% Compute observed pointwise F-statistics
    [fStats, pValues, dfWithinAt, etaSqPartial] = compute_pointwise_F( ...
        cycles, blockLabel, blockIDs, nBlocks, minValidPerBlk, cycleLength);

    dfBetween = nBlocks - 1;

    %% Multiple comparisons correction
    pCorrected         = NaN(1, cycleLength);
    significant        = false(1, cycleLength);
    sigClusters        = struct('startIdx',{},'endIdx',{},'startTime',{}, ...
                                'endTime',{},'duration',{},'meanF',{}, ...
                                'peakF',{},'mass',{},'pValue',{});
    allClusters        = sigClusters;
    clusterFThreshold  = NaN;
    nullMaxClusterMass = [];

    switch correction
        case {'none', 'bonferroni', 'fdr'}
            %% Pointwise correction
            pCorrected = apply_pointwise_correction(pValues, correction);
            significant = pCorrected < alpha;
            significant(isnan(pCorrected)) = false;

            % Identify clusters for reporting (no permutation p-values)
            allClusters = find_and_label_clusters(significant, fStats, ...
                                                   timeAxis);
            sigClusters = allClusters;
            for c = 1:length(sigClusters)
                sigClusters(c).pValue = NaN;
            end

        case 'cluster'
            %% Cluster-based permutation testing
            % Determine cluster-forming F threshold
            medianDfW = median(dfWithinAt, 'omitnan');
            if isnan(medianDfW) || medianDfW <= 0
                warning('block_ftest:noValidDf', ...
                    'Cannot compute cluster threshold — insufficient data.');
                clusterFThreshold = Inf;
            else
                clusterFThreshold = finv(1 - clusterAlpha, dfBetween, ...
                                          medianDfW);
            end

            % Form observed clusters
            obsSupra = fStats > clusterFThreshold & ~isnan(fStats);
            allClusters = find_and_label_clusters(obsSupra, fStats, ...
                                                   timeAxis);

            if isempty(allClusters)
                nullMaxClusterMass = [];
            else
                % Permutation null distribution
                fprintf('  block_ftest: Running %d permutations', ...
                    nPermutations);
                nullMaxClusterMass = zeros(nPermutations, 1);
                reportInterval = max(round(nPermutations / 10), 1);

                for perm = 1:nPermutations
                    if mod(perm, reportInterval) == 0
                        fprintf('.');
                    end

                    % Shuffle block labels
                    shuffIdx = randperm(nCycles);
                    shuffLabels = blockLabel(shuffIdx);

                    % Recompute F-statistics (no p-values needed)
                    permF = compute_pointwise_F_fast(cycles, ...
                                shuffLabels, blockIDs, nBlocks, ...
                                minValidPerBlk, cycleLength);

                    % Find largest cluster mass under permutation
                    permSupra = permF > clusterFThreshold & ~isnan(permF);
                    if any(permSupra)
                        permClusters = find_and_label_clusters( ...
                                           permSupra, permF, timeAxis);
                        if ~isempty(permClusters)
                            nullMaxClusterMass(perm) = ...
                                max([permClusters.mass]);
                        end
                    end
                end
                fprintf(' done.\n');

                % Assign cluster-level p-values
                for c = 1:length(allClusters)
                    allClusters(c).pValue = mean( ...
                        nullMaxClusterMass >= allClusters(c).mass);
                end

                % Identify significant clusters
                sigMask = [allClusters.pValue] < alpha;
                sigClusters = allClusters(sigMask);

                % Build pointwise significance mask
                for c = 1:length(sigClusters)
                    idx = sigClusters(c).startIdx:sigClusters(c).endIdx;
                    significant(idx) = true;
                end
            end
    end

    %% Summary statistics
    nSigPoints   = sum(significant);
    nValidTests  = sum(~isnan(pValues));
    fracSigPoints = nSigPoints / max(nValidTests, 1);

    [maxF, maxFIdx] = max(fStats);
    if ~isempty(maxFIdx) && ~isnan(maxF)
        maxFTime = timeAxis(maxFIdx);
    else
        maxFTime = NaN;
    end

    recommendation = generate_recommendation(nBlocks, nSigPoints, ...
        nValidTests, fracSigPoints, sigClusters, allClusters, ...
        blockIDs, correction);

    %% Post-hoc pairwise comparisons
    pairs  = nchoosek(1:nBlocks, 2);
    nPairs = size(pairs, 1);

    pairTStats    = NaN(nPairs, cycleLength);
    pairPValues   = NaN(nPairs, cycleLength);
    pairPCorrPost = NaN(nPairs, cycleLength);
    pairSig       = false(nPairs, cycleLength);
    pairMeanDiffs = NaN(nPairs, cycleLength);

    sigIdx = find(significant);

    for t = sigIdx
        for pp = 1:nPairs
            b1 = blockIDs(pairs(pp, 1));
            b2 = blockIDs(pairs(pp, 2));

            y1 = cycles(blockLabel == b1, t);
            y2 = cycles(blockLabel == b2, t);
            y1 = y1(~isnan(y1));
            y2 = y2(~isnan(y2));

            if length(y1) < 2 || length(y2) < 2
                continue;
            end

            n1 = length(y1);
            n2 = length(y2);
            m1 = mean(y1);
            m2 = mean(y2);

            pairMeanDiffs(pp, t) = m1 - m2;

            s1 = var(y1);
            s2 = var(y2);
            se = sqrt(s1/n1 + s2/n2);

            if se > 0
                tStat = (m1 - m2) / se;
                num = (s1/n1 + s2/n2)^2;
                den = (s1/n1)^2/(n1-1) + (s2/n2)^2/(n2-1);
                df  = num / den;

                pairTStats(pp, t)  = tStat;
                pairPValues(pp, t) = 2 * (1 - tcdf(abs(tStat), df));
            end
        end

        rawP = pairPValues(:, t);
        switch postHocMethod
            case 'holm'
                pairPCorrPost(:, t) = holm_bonferroni(rawP);
            case 'tukey'
                balanced = all(nCyclesPerBlock == nCyclesPerBlock(1));
                if balanced
                    pairPCorrPost(:, t) = tukey_correction( ...
                        rawP, nPairs, nBlocks, ...
                        median(dfWithinAt, 'omitnan'));
                else
                    pairPCorrPost(:, t) = holm_bonferroni(rawP);
                end
        end

        pairSig(:, t) = pairPCorrPost(:, t) < alpha;
    end

    pairBlockIDs = zeros(nPairs, 2);
    for pp = 1:nPairs
        pairBlockIDs(pp, :) = [blockIDs(pairs(pp,1)), ...
                                blockIDs(pairs(pp,2))];
    end

    pairwise = struct();
    pairwise.pairs       = pairBlockIDs;
    pairwise.tStats      = pairTStats;
    pairwise.pValues     = pairPValues;
    pairwise.pCorrected  = pairPCorrPost;
    pairwise.significant = pairSig;
    pairwise.meanDiffs   = pairMeanDiffs;

    %% Pack output
    results = struct();

    % --- Metadata ---
    results.nBlocks            = nBlocks;
    results.blockIDs           = blockIDs;
    results.nCyclesPerBlock    = nCyclesPerBlock;
    results.cycleLength        = cycleLength;
    results.alpha              = alpha;
    results.correction         = correction;
    results.postHocMethod      = postHocMethod;
    results.timeAxis           = timeAxis;

    % --- Block-level descriptives ---
    results.blockMeans         = blockMeans;
    results.blockSEMs          = blockSEMs;
    results.blockNValid        = blockNValid;
    results.grandMean          = grandMean;

    % --- Omnibus F-test ---
    results.fStats             = fStats;
    results.pValues            = pValues;
    results.pCorrected         = pCorrected;
    results.significant        = significant;
    results.dfBetween          = dfBetween;
    results.dfWithinAt         = dfWithinAt;
    results.etaSquaredPartial  = etaSqPartial;

    % --- Summary ---
    results.nSigPoints         = nSigPoints;
    results.fracSigPoints      = fracSigPoints;
    results.maxF               = maxF;
    results.maxFTime           = maxFTime;
    results.sigClusters        = sigClusters;
    results.allClusters        = allClusters;
    results.recommendation     = recommendation;

    % --- Cluster permutation details ---
    results.clusterAlpha       = clusterAlpha;
    results.clusterFThreshold  = clusterFThreshold;
    results.nPermutations      = nPermutations;
    results.nullMaxClusterMass = nullMaxClusterMass;

    % --- Post-hoc ---
    results.pairwise           = pairwise;

end


%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%  ========================================================================

function [fStats, pValues, dfWithinAt, etaSqPartial] = ...
    compute_pointwise_F(cycles, blockLabel, blockIDs, nBlocks, ...
                        minValidPerBlk, cycleLength)
    % COMPUTE_POINTWISE_F  Vectorized one-way ANOVA F-statistic at each time point.

    dfB = nBlocks - 1;

    % Build a [nCycles x nBlocks] indicator matrix
    % membership(i,b) = 1 if cycle i belongs to block b
    membership = blockLabel == blockIDs';  % [nCycles x nBlocks]

    % Valid (non-NaN) mask: [nCycles x cycleLength]
    validMask = ~isnan(cycles);

    % Replace NaN with 0 for summation (NaN would propagate)
    cyclesZeroed = cycles;
    cyclesZeroed(~validMask) = 0;

    % Per-block valid counts at each timepoint: [nBlocks x cycleLength]
    %   blockN(b,t) = number of valid cycles in block b at timepoint t
    blockN = membership' * validMask;  % matrix multiply: [nBlocks x nCycles] * [nCycles x cycleLength]

    % Per-block sums at each timepoint: [nBlocks x cycleLength]
    blockSum = membership' * (cyclesZeroed .* validMask);

    % Per-block means: [nBlocks x cycleLength]
    blockMean = blockSum ./ max(blockN, 1);

    % Total valid count and grand sum per timepoint: [1 x cycleLength]
    nTotal = sum(blockN, 1);
    grandSum = sum(blockSum, 1);
    grandMean = grandSum ./ max(nTotal, 1);

    % --- SS Between: sum_b [ n_b * (blockMean_b - grandMean)^2 ] ---
    ssBetween = sum(blockN .* (blockMean - grandMean).^2, 1);  % [1 x cycleLength]

    % --- SS Within: sum of (x_i - blockMean_b)^2 for each cycle i in block b ---
    %   Trick: SS_within = sum(x^2) - sum_b[ blockSum_b^2 / n_b ]
    %   This avoids expanding blockMean back to per-cycle resolution.
    sumSq = membership' * (cyclesZeroed.^2 .* validMask);  % [nBlocks x cycleLength]
    ssWithin = sum(sumSq - blockSum.^2 ./ max(blockN, 1), 1);  % [1 x cycleLength]

    % --- Degrees of freedom ---
    nBlocksPresent = sum(blockN >= 1, 1);  % [1 x cycleLength]
    dfW = nTotal - nBlocksPresent;         % [1 x cycleLength]

    % --- Validity mask: all blocks must meet minValidPerBlock ---
    allBlocksValid = all(blockN >= minValidPerBlk, 1);  % [1 x cycleLength]
    computable = allBlocksValid & nBlocksPresent >= 2 & dfW > 0;

    % --- Compute F-statistics ---
    msBetween = ssBetween / dfB;
    msWithin  = ssWithin ./ max(dfW, 1);

    fStats       = NaN(1, cycleLength);
    pValues      = NaN(1, cycleLength);
    dfWithinAt   = NaN(1, cycleLength);
    etaSqPartial = NaN(1, cycleLength);

    validF = computable & msWithin > 0;
    fStats(validF) = msBetween(validF) ./ msWithin(validF);
    dfWithinAt(computable) = dfW(computable);

    % p-values (fcdf is vectorized)
    pValues(validF) = 1 - fcdf(fStats(validF), dfB, dfW(validF));

    % Handle zero within-variance
    zeroVar = computable & msWithin == 0;
    fStats(zeroVar) = Inf;
    pValues(zeroVar) = 0;

    % Effect size
    ssTotal = ssBetween + ssWithin;
    validEta = computable & ssTotal > 0;
    etaSqPartial(validEta) = ssBetween(validEta) ./ ssTotal(validEta);
end


function fStats = compute_pointwise_F_fast(cycles, blockLabel, ...
    blockIDs, nBlocks, minValidPerBlk, cycleLength)
    % COMPUTE_POINTWISE_F_FAST  Vectorized F-statistics only (for permutation speed).
    %   Stripped version that skips p-values, df, and effect size.

    dfB = nBlocks - 1;

    membership = blockLabel == blockIDs';
    validMask = ~isnan(cycles);
    cyclesZeroed = cycles;
    cyclesZeroed(~validMask) = 0;

    blockN   = membership' * validMask;
    blockSum = membership' * (cyclesZeroed .* validMask);
    blockMean = blockSum ./ max(blockN, 1);

    nTotal    = sum(blockN, 1);
    grandMean = sum(blockSum, 1) ./ max(nTotal, 1);

    ssBetween = sum(blockN .* (blockMean - grandMean).^2, 1);

    sumSq    = membership' * (cyclesZeroed.^2 .* validMask);
    ssWithin = sum(sumSq - blockSum.^2 ./ max(blockN, 1), 1);

    nBlocksPresent = sum(blockN >= 1, 1);
    dfW = nTotal - nBlocksPresent;

    allBlocksValid = all(blockN >= minValidPerBlk, 1);
    computable = allBlocksValid & nBlocksPresent >= 2 & dfW > 0;

    fStats = NaN(1, cycleLength);
    msWithin = ssWithin ./ max(dfW, 1);
    validF = computable & msWithin > 0;
    fStats(validF) = (ssBetween(validF) / dfB) ./ msWithin(validF);
end


function clusters = find_and_label_clusters(supraMask, fStats, timeAxis)
% FIND_AND_LABEL_CLUSTERS  Find contiguous supra-threshold regions.
    clusters = struct('startIdx',{},'endIdx',{},'startTime',{}, ...
                      'endTime',{},'duration',{},'meanF',{}, ...
                      'peakF',{},'mass',{},'pValue',{});

    if ~any(supraMask)
        return;
    end

    dt = 0;
    if length(timeAxis) > 1
        dt = timeAxis(2) - timeAxis(1);
    end

    d = diff([0, supraMask, 0]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;

    for c = 1:length(starts)
        idx = starts(c):ends(c);
        cl = struct();
        cl.startIdx  = starts(c);
        cl.endIdx    = ends(c);
        cl.startTime = timeAxis(starts(c));
        cl.endTime   = timeAxis(ends(c));
        cl.duration  = cl.endTime - cl.startTime + dt;
        cl.meanF     = mean(fStats(idx), 'omitnan');
        cl.peakF     = max(fStats(idx));
        cl.mass      = sum(fStats(idx), 'omitnan');
        cl.pValue    = NaN;
        clusters(end+1) = cl; %#ok<AGROW>
    end
end


function pCorr = apply_pointwise_correction(pRaw, method)
% APPLY_POINTWISE_CORRECTION  Bonferroni or FDR correction.
    pCorr = pRaw;
    validIdx = find(~isnan(pRaw));
    if isempty(validIdx)
        return;
    end
    pValid = pRaw(validIdx);

    switch method
        case 'none'
            % No correction
        case 'bonferroni'
            nTests = length(pValid);
            pValid = min(pValid * nTests, 1);
        case 'fdr'
            pValid = benjamini_hochberg(pValid);
    end

    pCorr(validIdx) = pValid;
end


function pFDR = benjamini_hochberg(pRaw)
% BENJAMINI_HOCHBERG  Benjamini-Hochberg FDR correction.
    m = length(pRaw);
    [pSorted, sortIdx] = sort(pRaw);
    pAdj = zeros(m, 1);

    for k = 1:m
        pAdj(k) = pSorted(k) * m / k;
    end

    for k = m-1:-1:1
        pAdj(k) = min(pAdj(k), pAdj(k+1));
    end

    pAdj = min(pAdj, 1);

    pFDR = zeros(m, 1);
    pFDR(sortIdx) = pAdj;
    pFDR = pFDR';
end


function pCorr = holm_bonferroni(pRaw)
% HOLM_BONFERRONI  Holm-Bonferroni step-down correction.
    pRaw = pRaw(:);
    m = length(pRaw);
    validIdx = find(~isnan(pRaw));
    pCorr = pRaw;

    if isempty(validIdx)
        return;
    end

    pValid = pRaw(validIdx);
    nValid = length(pValid);
    [pSorted, sortIdx] = sort(pValid);

    pAdj = zeros(nValid, 1);
    for k = 1:nValid
        pAdj(k) = pSorted(k) * (nValid - k + 1);
    end

    for k = 2:nValid
        pAdj(k) = max(pAdj(k), pAdj(k-1));
    end

    pAdj = min(pAdj, 1);

    pOut = zeros(nValid, 1);
    pOut(sortIdx) = pAdj;
    pCorr(validIdx) = pOut;
end


function pCorr = tukey_correction(pRaw, ~, nGroups, ~)
% TUKEY_CORRECTION  Approximate Tukey HSD correction (Sidak-like).
    pRaw = pRaw(:);
    pCorr = pRaw;
    nComparisons = nchoosek(nGroups, 2);

    for k = 1:length(pRaw)
        if isnan(pRaw(k))
            continue;
        end
        pCorr(k) = min(1 - (1 - pRaw(k))^nComparisons, 1);
    end

    pCorr = min(max(pCorr, 0), 1);
end


function txt = generate_recommendation(nBlocks, nSig, nValid, fracSig, ...
                                        sigClusters, allClusters, ...
                                        blockIDs, correction)
% GENERATE_RECOMMENDATION  Create summary text from F-test results.
    nSigC = length(sigClusters);
    nAllC = length(allClusters);

    if strcmp(correction, 'cluster')
        corrStr = sprintf('cluster permutation, %d/%d clusters significant', ...
            nSigC, nAllC);
    else
        corrStr = sprintf('%s corrected', correction);
    end

    if nSig == 0
        txt = sprintf(['No significant block effects detected (%s). ' ...
            'The %d blocks appear exchangeable — pooling all cycles ' ...
            'is appropriate.'], corrStr, nBlocks);
    elseif fracSig < 0.05
        txt = sprintf(['Minor block effects at %.1f%% of time points ' ...
            '(%d/%d, %s). Blocks are largely exchangeable; pooling is ' ...
            'likely acceptable, but inspect the affected regions.'], ...
            fracSig*100, nSig, nValid, corrStr);
    elseif fracSig < 0.30
        if nSigC > 0
            durStr = sprintf('%.0f-%.0f ms', ...
                min([sigClusters.duration])*1000, ...
                max([sigClusters.duration])*1000);
        else
            durStr = 'N/A';
        end
        txt = sprintf(['Moderate block effects at %.1f%% of time points ' ...
            '(%d/%d, %s). %d significant cluster(s) spanning %s. ' ...
            'Inspect block means to determine if effects are ' ...
            'scientifically meaningful or confounding.'], ...
            fracSig*100, nSig, nValid, corrStr, nSigC, durStr);
    else
        txt = sprintf(['Substantial block effects at %.1f%% of time ' ...
            'points (%d/%d, %s). Blocks are NOT exchangeable. ' ...
            'Do not pool — average block means instead, or model the ' ...
            'block effect explicitly.'], ...
            fracSig*100, nSig, nValid, corrStr);
    end
end