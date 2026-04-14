function timepoints = run_timepointAnalysis(blocks, params, cycleLength, nPermutations)
% RUN_TIMEPOINTANALYSIS  Groups blocks by timepoint label and computes
% pooled cycle averages, cycle-metrics, and the Omnibus F-test.

    if nargin < 4, nPermutations = params.ftestPerturbs; end
    timepoints = struct;
    tp_labels = string({blocks.timePoint});
    tp_labels = tp_labels(~strcmp(tp_labels,'NA'));
    tp_unique = unique(tp_labels);
    ntp = length(tp_unique);

    for ii = 1:ntp
        tp_ii = tp_unique(ii);
        tp_ii_ids = find(strcmp({blocks.timePoint}, tp_ii));
        tp_blocknumbers = strjoin(string([blocks(tp_ii_ids).blockNumber]), ', ');
        nids = length(tp_ii_ids);

        timepoint_cyclelabels = [];
        good_cyclemean_mat   = zeros(nids, cycleLength);
        good_cyclemedian_mat = zeros(nids, cycleLength);

        for jj = 1:nids
            id_jj = tp_ii_ids(jj);
            ngc_jj = size(blocks(id_jj).good_cyclemat, 1);
            timepoint_cyclelabels = [timepoint_cyclelabels; repmat(id_jj, ngc_jj, 1)];
            good_cyclemean_mat(jj,:)   = mean(blocks(id_jj).good_cyclemat, 1, 'omitnan');
            good_cyclemedian_mat(jj,:) = median(blocks(id_jj).good_cyclemat, 1, 'omitnan');
        end

        drumvel_cyclemat   = vertcat(blocks(tp_ii_ids).stimvel_cyclemat);
        timepoint_cyclemat = vertcat(blocks(tp_ii_ids).good_cyclemat);

        res_ii = block_ftest(timepoint_cyclemat, timepoint_cyclelabels, ...
            'Fs',params.fs, 'StimFreq',params.exp_stimfreq, ...
            'nPermutations',nPermutations);

        nPooledCycles = size(timepoint_cyclemat, 1);
        pooled_drumvel_cyclemean_cvd = cycleVarianceDecomposition(drumvel_cyclemat);
        pooled_good_cyclemean_cvd    = cycleVarianceDecomposition(timepoint_cyclemat);
        pooled_good_cyclemedian      = median(timepoint_cyclemat, 1, 'omitnan');

        good_cyclemean_cvd   = cycleVarianceDecomposition(good_cyclemean_mat);
        good_cyclemedian_cvd = cycleVarianceDecomposition(good_cyclemedian_mat);

        nGoodCycles = [blocks(tp_ii_ids).nGoodCycles];
        good_nGoodCyclesWeighted_cyclemean = ...
            sum(good_cyclemean_mat.*nGoodCycles', 1) / sum(nGoodCycles');

        check_valuesApproxEqual(pooled_good_cyclemean_cvd.cycleMean, ...
            good_nGoodCyclesWeighted_cyclemean, ...
            [char(tp_ii), ' cycle-means']);

        timepoints(ii).timePoint        = char(tp_ii);
        timepoints(ii).blockType        = blocks(tp_ii_ids(1)).blockType;
        timepoints(ii).blockNumbers     = tp_blocknumbers;
        timepoints(ii).blockResultsIds  = tp_ii_ids;

        pooled_good_cyclemean_cm        = calc_cycleMetrics(pooled_good_cyclemean_cvd.cycleMean, pooled_drumvel_cyclemean_cvd.cycleMean);
        pooled_good_cyclemedian_cm      = calc_cycleMetrics(pooled_good_cyclemedian,              pooled_drumvel_cyclemean_cvd.cycleMean);
        good_cyclemean_cm               = calc_cycleMetrics(good_cyclemean_cvd.cycleMean,         pooled_drumvel_cyclemean_cvd.cycleMean);
        good_cyclemedian_cm             = calc_cycleMetrics(good_cyclemedian_cvd.cycleMean,       pooled_drumvel_cyclemean_cvd.cycleMean);
        nGoodCyclesWeighted_cyclemean_cm = calc_cycleMetrics(good_nGoodCyclesWeighted_cyclemean,  pooled_drumvel_cyclemean_cvd.cycleMean);

        timepoints(ii).nPooledCycles                    = nPooledCycles;
        timepoints(ii).pooled_drumvel_cyclemean_cvd     = pooled_drumvel_cyclemean_cvd;
        timepoints(ii).pooled_good_cyclemean_cvd        = pooled_good_cyclemean_cvd;
        timepoints(ii).pooled_good_cyclemean_cm         = pooled_good_cyclemean_cm;
        timepoints(ii).pooled_good_cyclemedian          = pooled_good_cyclemedian;
        timepoints(ii).pooled_good_cyclemedian_cm       = pooled_good_cyclemedian_cm;

        timepoints(ii).good_cyclemean_cvd               = good_cyclemean_cvd;
        timepoints(ii).good_cyclemean_cm                = good_cyclemean_cm;
        timepoints(ii).good_cyclemedian_cvd             = good_cyclemedian_cvd;
        timepoints(ii).good_cyclemedian_cm              = good_cyclemedian_cm;
        timepoints(ii).nGoodCyclesWeighted_cyclemean    = good_nGoodCyclesWeighted_cyclemean;
        timepoints(ii).nGoodCyclesWeighted_cyclemean_cm = nGoodCyclesWeighted_cyclemean_cm;

        timepoints(ii).good_cyclemean_mat   = good_cyclemean_mat;
        timepoints(ii).good_cyclemedian_mat = good_cyclemedian_mat;

        timepoints(ii).ftestFracSigPoints = res_ii.fracSigPoints;
        timepoints(ii).ftestEtaSquared    = res_ii.etaSquaredPartial;
        timepoints(ii).ftestStats         = res_ii.fStats;
    end
end

function check_valuesApproxEqual(a, b, label, tol)
    if ~exist('tol', 'var'); tol = 5e-7; end
    if any(isnan([a(:); b(:)]))
        warning('NaN values for %s cannot be compared', label);
        return;
    end
    if ~all(abs(a-b) < tol)
        warning('Values for %s do not agree', label);
    end
end
