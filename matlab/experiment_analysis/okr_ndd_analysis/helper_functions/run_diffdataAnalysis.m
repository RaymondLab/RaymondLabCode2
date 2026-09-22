function diffdata = run_diffdataAnalysis(timepoints)
% RUN_DIFFDATAANALYSIS  Computes post-minus-pre cycle-mean differences and
% cycle metrics for every ordered pair of timepoints.

    diffdata = struct;
    ntp = length(timepoints);
    tp_pair_ids = get_diffdata_permutations(1:ntp);
    npairs = size(tp_pair_ids, 1);

    for ii = 1:npairs
        preidx  = tp_pair_ids(ii,1);
        postidx = tp_pair_ids(ii,2);
        diffdata(ii).timePointDiff        = sprintf('%s-%s', timepoints(postidx).timePoint, timepoints(preidx).timePoint);
        diffdata(ii).timePointPrePostIds  = [preidx, postidx];

        diffdata(ii).good_cyclemean_diff = ...
            timepoints(postidx).good_cyclemean_cvd.cycleMean - timepoints(preidx).good_cyclemean_cvd.cycleMean;
        diffdata(ii).good_cyclemean_diff_SEM = ...
            calc_cyclemeanDiffSEM(timepoints(preidx).good_cyclemean_cvd.cycleSEM, timepoints(postidx).good_cyclemean_cvd.cycleSEM);
        diffdata(ii).good_cyclemean_diff_cm = ...
            calc_cycleMetrics(diffdata(ii).good_cyclemean_diff, timepoints(preidx).pooled_drumvel_cyclemean_cvd.cycleMean);

        diffdata(ii).nGoodCyclesWeighted_cyclemean_diff = ...
            timepoints(postidx).nGoodCyclesWeighted_cyclemean - timepoints(preidx).nGoodCyclesWeighted_cyclemean;
        diffdata(ii).nGoodCyclesWeighted_cyclemean_diff_cm = ...
            calc_cycleMetrics(diffdata(ii).nGoodCyclesWeighted_cyclemean_diff, timepoints(preidx).pooled_drumvel_cyclemean_cvd.cycleMean);

        diffdata(ii).good_cyclemedianmean_diff = ...
            timepoints(postidx).good_cyclemedian_cvd.cycleMean - timepoints(preidx).good_cyclemedian_cvd.cycleMean;
        diffdata(ii).good_cyclemedianmean_diff_cm = ...
            calc_cycleMetrics(diffdata(ii).good_cyclemedianmean_diff, timepoints(preidx).pooled_drumvel_cyclemean_cvd.cycleMean);

        diffdata(ii).pooled_good_cyclemean_diff = ...
            timepoints(postidx).pooled_good_cyclemean_cvd.cycleMean - timepoints(preidx).pooled_good_cyclemean_cvd.cycleMean;
        diffdata(ii).pooled_good_cyclemean_diff_cm = ...
            calc_cycleMetrics(diffdata(ii).pooled_good_cyclemean_diff, timepoints(preidx).pooled_drumvel_cyclemean_cvd.cycleMean);
        diffdata(ii).pooled_good_cyclemean_diff_SEM = ...
            calc_cyclemeanDiffSEM(timepoints(preidx).pooled_good_cyclemean_cvd.cycleSEM, timepoints(postidx).pooled_good_cyclemean_cvd.cycleSEM);

        diffdata(ii).pooled_good_cyclemedian_diff = ...
            timepoints(postidx).pooled_good_cyclemedian - timepoints(preidx).pooled_good_cyclemedian;
        diffdata(ii).pooled_good_cyclemedian_diff_cm = ...
            calc_cycleMetrics(diffdata(ii).pooled_good_cyclemedian_diff, timepoints(preidx).pooled_drumvel_cyclemean_cvd.cycleMean);
    end
end

function diffSEM = calc_cyclemeanDiffSEM(cyclesem1, cyclesem2)
    diffSEM = sqrt(cyclesem1.^2 + cyclesem2.^2);
end
