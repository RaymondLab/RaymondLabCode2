function block = run_blockAnalysis(params, data, ii, lpc, saccThresh)
% RUN_BLOCKANALYSIS  Runs the per-block OKR analysis for block index ii,
% parameterized on lowpass cutoff (lpc) and saccade detection threshold
% (saccThresh). Returns a single struct row (one element of `blocks`).
% The caller is responsible for setting block.timePoint.

    blockIds    = data.blockIds;
    blockTypes  = data.blockTypes;
    bssTimes    = data.bssTimes;
    cycleLength = data.cycleLength;
    cycleTimes  = data.cycleTimes;
    chairvel_raw = data.chairvel_raw;
    drumvel_raw  = data.drumvel_raw;
    hevel_raw    = data.hevel_raw;
    hepos_raw    = data.hepos_raw;
    saccadePad   = data.saccadePad;

    % Extract corresponding block segments for chair, drum, and eye data
    chairvel_raw_ii = chairvel_raw(blockIds(1,ii):blockIds(2,ii));
    drumvel_raw_ii  = drumvel_raw(blockIds(1,ii):blockIds(2,ii));
    hevel_raw_ii    = hevel_raw(blockIds(1,ii):blockIds(2,ii));
    hepos_raw_ii    = hepos_raw(blockIds(1,ii):blockIds(2,ii));

    blockLength = length(drumvel_raw_ii);
    blockTimes  = (0:blockLength-1) / params.fs;

    % Perform initial fit of the raw eye velocity data
    keep = abs(hevel_raw_ii) < 5*std(abs(hevel_raw_ii)) + mean(abs(hevel_raw_ii));
    [fit0,~] = calc_sineFit(params.exp_stimfreq, params.fs, keep, hevel_raw_ii);
    hevel_raw_fit = fit0.eyevel_fit;

    % Lowpass-filter eye position, take moving-slope to get velocity
    hepos_ii = butterworthfilter(hepos_raw_ii, lpc, params.fs, 'n',4);
    hevel_ii = movingslope(hepos_ii, params.filterWindow) * params.fs;

    % Apply saccade detection
    [saccDilMask, saccMask, hevel_des_ii] = run_saccadeDetection(hevel_ii, ...
        hevel_raw_fit, saccThresh, saccadePad, ...
        params.minGoodChunk_len, params.saccadeMethod);

    [bfit,stat] = calc_sineFit(params.exp_stimfreq, params.fs, ~saccDilMask, ...
        hevel_ii, chairvel_raw_ii, drumvel_raw_ii);

    eyevel_des_cyclefit = sin(2*pi*params.exp_stimfreq*cycleTimes + ...
        deg2rad(bfit.eyevel_rel_phase+180)) * bfit.eyevel_amp;

    startpt = max(1, round(mod(-bfit.ref_angle,360) / 360 * params.fs/params.exp_stimfreq));

    cyclestart_ids_mat = segmentIntoCycles((1:blockLength), startpt, cycleLength);
    startpts = cyclestart_ids_mat(:,1);

    chairvel_cyclemat  = segmentIntoCycles(chairvel_raw_ii, startpt, cycleLength);
    drumvel_cyclemat   = segmentIntoCycles(drumvel_raw_ii, startpt, cycleLength);
    hevel_cyclemat     = segmentIntoCycles(hevel_ii, startpt, cycleLength);
    hevel_des_cyclemat = segmentIntoCycles(hevel_des_ii, startpt, cycleLength);
    saccmask_mat       = segmentIntoCycles(double(saccMask), startpt, cycleLength);
    mse_cyclemat       = segmentIntoCycles((hevel_ii-hevel_raw_fit).^2, startpt, cycleLength);

    badCycles        = any(saccmask_mat, 2);
    nGoodCycles      = sum(~badCycles);                  % Number of "good" cycles (clean of saccades)
    [nTotalCycles,~] = size(hevel_cyclemat);             % Number of total cycles in the block
    saccadeFrac      = sum(saccMask) / numel(saccMask);  % Fraction of block samples that were flagged as saccades

    chairvel_cyclecvd = cycleVarianceDecomposition(chairvel_cyclemat);
    drumvel_cyclecvd  = cycleVarianceDecomposition(drumvel_cyclemat);
    hevel_cyclecvd    = cycleVarianceDecomposition(hevel_cyclemat);

    if nGoodCycles > 0
        hevel_good_cyclemat = hevel_cyclemat(~badCycles, :);
        hevel_good_cyclecvd = cycleVarianceDecomposition(hevel_good_cyclemat);
    else
        hevel_good_cyclemat = zeros(size(hevel_cyclemat));
        hevel_good_cyclecvd = cycleVarianceDecomposition(zeros(size(hevel_cyclemat(1,:))));
    end

    [cfit,~] = calc_sineFit(params.exp_stimfreq, params.fs, [], ...
        hevel_good_cyclecvd.cycleMean, chairvel_cyclecvd.cycleMean, drumvel_cyclecvd.cycleMean);

    % Populate block result struct
    block = struct;
    block.blockType   = char(blockTypes(ii));
    block.blockNumber = ii;

    if bfit.chairvel_amp > 3
        block.stimType = 'Chair';
        [good_rel_gain,good_rel_gain_SEM] = calc_eyeRelGainSEM(hevel_good_cyclecvd.amplitudeMean, ...
            hevel_good_cyclecvd.amplitudeSEM, ...
            chairvel_cyclecvd.amplitudeMean, chairvel_cyclecvd.amplitudeSEM);
        [good_rel_phase,good_rel_phase_SEM] = calc_eyeRelPhaseSEM(hevel_good_cyclecvd.phaseMean_deg, ...
            hevel_good_cyclecvd.phaseSEM_deg, ...
            chairvel_cyclecvd.phaseMean_deg, chairvel_cyclecvd.phaseSEM_deg);
    elseif bfit.drumvel_amp > 3
        block.stimType = 'Drum';
        [good_rel_gain,good_rel_gain_SEM] = calc_eyeRelGainSEM(hevel_good_cyclecvd.amplitudeMean, ...
            hevel_good_cyclecvd.amplitudeSEM, ...
            drumvel_cyclecvd.amplitudeMean, drumvel_cyclecvd.amplitudeSEM);
        [good_rel_phase,good_rel_phase_SEM] = calc_eyeRelPhaseSEM(hevel_good_cyclecvd.phaseMean_deg, ...
            hevel_good_cyclecvd.phaseSEM_deg, ...
            drumvel_cyclecvd.phaseMean_deg, drumvel_cyclecvd.phaseSEM_deg);
    else
        block.stimType = 'NoStim';
        good_rel_gain = nan; good_rel_gain_SEM = nan;
        good_rel_phase = nan; good_rel_phase_SEM = nan;
    end

    block.stimFreq = params.exp_stimfreq;
    block.hevel_amp       = bfit.eyevel_amp;
    block.hevel_phase     = bfit.eyevel_phase;
    block.hevel_rel_gain  = bfit.eyevel_rel_gain;
    block.hevel_rel_phase = bfit.eyevel_rel_phase;
    block.good_amp           = hevel_good_cyclecvd.amplitudeMean;
    block.good_amp_SEM       = hevel_good_cyclecvd.amplitudeSEM;
    block.good_phase         = hevel_good_cyclecvd.phaseMean_deg;
    block.good_phase_SEM     = hevel_good_cyclecvd.phaseSEM_deg;
    block.good_rel_gain      = good_rel_gain;
    block.good_rel_gain_SEM  = good_rel_gain_SEM;
    block.good_rel_phase     = good_rel_phase;
    block.good_rel_phase_SEM = good_rel_phase_SEM;
    block.saccadeFrac  = saccadeFrac;
    block.nGoodCycles  = nGoodCycles;
    block.nTotalCycles = nTotalCycles;

    block.startTime = bssTimes(1,ii);
    block.endTime   = bssTimes(2,ii);
    if bfit.chairvel_amp > 3
        block.stimvel = chairvel_raw_ii;
    elseif bfit.drumvel_amp > 3
        block.stimvel = drumvel_raw_ii;
    else
        block.stimvel = zeros(size(blockTimes));
    end
    block.hepos_raw     = hepos_raw_ii;
    block.hevel_raw     = hevel_raw_ii;
    block.hevel         = hevel_ii;
    block.hevel_des     = hevel_des_ii;
    block.hevel_raw_fit = hevel_raw_fit;
    block.hevel_des_fit = bfit.eyevel_fit;

    block.cycleStartIds = startpts;
    block.des_cyclefit  = eyevel_des_cyclefit;
    block.good_cyclefit = cfit.eyevel_rel_fit;
    if bfit.chairvel_amp > 3
        block.stimvel_cyclemat = chairvel_cyclemat;
    elseif bfit.drumvel_amp > 3
        block.stimvel_cyclemat = drumvel_cyclemat;
    else
        block.stimvel_cyclemat = zeros(size(hevel_cyclemat));
    end
    block.cyclemat      = hevel_cyclemat;
    block.des_cyclemat  = hevel_des_cyclemat;
    block.good_cyclemat = hevel_good_cyclemat;
    block.mse_cyclemat  = mse_cyclemat;

    block.headamp   = bfit.chairvel_amp;
    block.headangle = bfit.chairvel_angle;
    block.drumamp   = bfit.drumvel_amp;
    block.drumangle = bfit.drumvel_angle;
    block.rsquare   = stat(1);
    block.pval      = stat(3);
    block.variance  = sqrt(stat(4));

    if bfit.chairvel_amp > 3
        block.stimvel_cyclecvd = chairvel_cyclecvd;
        stimvel_amp = cfit.chairvel_amp;
        stimvel_phase = cfit.chairvel_angle;
    elseif bfit.drumvel_amp > 3
        block.stimvel_cyclecvd = drumvel_cyclecvd;
        stimvel_amp = cfit.drumvel_amp;
        stimvel_phase = cfit.drumvel_angle;
    else
        block.stimvel_cyclecvd = struct;
        stimvel_amp = nan; stimvel_phase = nan;
    end
    if ~isempty(fieldnames(block.stimvel_cyclecvd))
        block.hevel_cyclecm      = calc_cycleMetrics(hevel_cyclecvd.cycleMean, block.stimvel_cyclecvd.cycleMean);
        block.hevel_good_cyclecm = calc_cycleMetrics(hevel_good_cyclecvd.cycleMean, block.stimvel_cyclecvd.cycleMean);
    end
    block.hevel_cyclecvd      = hevel_cyclecvd;
    block.hevel_good_cyclecvd = hevel_good_cyclecvd;
    block.stimvel_amp       = block.stimvel_cyclecvd.amplitudeMean;
    block.stimvel_amp_SEM   = block.stimvel_cyclecvd.amplitudeSEM;
    block.stimvel_phase     = block.stimvel_cyclecvd.phaseMean_deg;
    block.stimvel_phase_SEM = block.stimvel_cyclecvd.phaseSEM_deg;
    block.amplitudeCV       = block.hevel_good_cyclecvd.amplitudeCV;
    block.residualMAD       = block.hevel_good_cyclecvd.residualMAD;
    block.etaSquared        = block.hevel_good_cyclecvd.etaSquared;

    check_valuesApproxEqual(cfit.eyevel_amp, block.hevel_good_cyclecvd.amplitudeMean, '"Good" cycles amplitudes');
    check_valuesApproxEqual(cfit.eyevel_phase, block.hevel_good_cyclecvd.phaseMean_deg, '"Good" cycles phases');
    check_valuesApproxEqual(stimvel_amp, block.stimvel_cyclecvd.amplitudeMean, 'Stim cycles amplitudes');
    check_valuesApproxEqual(stimvel_phase, block.stimvel_cyclecvd.phaseMean_deg, 'Stim cycles phases');
    check_valuesApproxEqual(cfit.eyevel_rel_gain, good_rel_gain, 'Relative gain');
    check_valuesApproxEqual(cfit.eyevel_rel_phase, good_rel_phase, 'Relative phase');
end

function check_valuesApproxEqual(a, b, label, tol)
    if ~exist('tol', 'var'); tol = 5e-7; end
    if any(isnan([a,b]))
        warning('NaN values for %s cannot be compared', label);
        return;
    end
    if ~all(abs(a-b) < tol)
        warning('Values for %s do not agree', label);
    end
end
