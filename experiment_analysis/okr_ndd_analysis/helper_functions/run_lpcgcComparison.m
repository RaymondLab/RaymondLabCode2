function [T, D] = run_lpcgcComparison(params, data)
% RUN_LPCGCCOMPARISON  Sweeps the analysis over the Cartesian product of
% params.lowpass_cutoffs and params.saccade_threshs. Returns cell arrays
% T{ii,jj} (timepoints) and D{ii,jj} (diffdata) for each grid point.

    ncutoff = length(params.lowpass_cutoffs);
    nthresh = length(params.saccade_threshs);
    T = cell(ncutoff, nthresh);
    D = cell(ncutoff, nthresh);

    nBlocks = length(params.timepoint_ids);

    for ii = 1:ncutoff
        lpc = params.lowpass_cutoffs(ii);
        for jj = 1:nthresh
            saccThresh = params.saccade_threshs(jj);
            fprintf('        LPC = %g Hz, Saccade Thresh = %g\n', lpc, saccThresh);

            clear blocks;
            for nn = 1:nBlocks
                blkIdx = params.timepoint_ids(nn);
                bgrp   = params.timepoint_groups(nn);
                bk = run_blockAnalysis(params, data, blkIdx, lpc, saccThresh);
                bk.timePoint = sprintf('T%d', bgrp);
                blocks(nn) = bk;
            end

            tpjj = run_timepointAnalysis(blocks, params, data.cycleLength, 100);
            ddjj = run_diffdataAnalysis(tpjj);
            T{ii,jj} = tpjj;
            D{ii,jj} = ddjj;
        end
    end
end
