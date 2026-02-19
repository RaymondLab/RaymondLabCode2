clear; clc; close all;
params.script_folderpath = fileparts(which(matlab.desktop.editor.getActiveFilename));
addpath(genpath(params.script_folderpath));

% OKR_NDD_EXP_ANALYSIS - OKR Normal-Dim-Dimmer (NDD) Experiment Analysis
% Important: This script requires "RaymondLabCode2" and its subfolders added to MATLAB path!


%% SCRIPT PARAMETERS (Do not edit unless you know what you are doing!)

% General experiment parameters
params.exp_stimfreq  = 1.0;    % Training stimulus frequency (Default 1.0)
params.exp_traintype = 'OKR';  % Training type (Default 'OKR')

% Define timepoint blocks and corresponding group ids
% Do not adjust this for now as the group_NDD_okrAnalysis is not set up for this
params.timepoint_ids    = [2, 3, 4, 59, 60, 61];  % Default [2, 3, 4, 59, 60, 61]
params.timepoint_groups = [0, 0, 0, 60, 60, 60];  % Default [0, 0, 0, 60, 60, 60]

% Subject condition options for dropdown menu
params.subcond_options = {'None', 'NR', 'DR'}; 

% Task (training) condition options for dropdown menu
params.taskcond_options = {'None', 'TNOR', 'TDIM', 'TDMR', '2TNOR', '2TDIM', '2TDMR'}; 

% Number of perturbations to run for F-Tests
params.ftestPerturbs = 10000;  % (Default 10000, or 5000 for tests)

% Custom channel scaling is used if tuple is provided (Default [])
params.customScaleChs = [];

% Data preprocessing parameters
params.filtering_options = {'Butter', 'FIF'};
params.transientThresh   = 1;   % Transient "spike" removal threshold (Default 1)
params.lowpassCutoff     = 9;   % Lowpass filter cutoff frequency (Default 11)
params.filterWindow      = 11;  % SG filter window size (Default 30)

% FIF-specific settings
params.fifSettings.delta = 0.1;  
params.fifSettings.Xi = 2;  
params.fifSettings.NIMFs = 4;  
params.fifSettings.alpha = 24;
params.fifSettings.minIMF = 4;

% Saccade removal parameters
params.saccadeMethod    = 'SVT';  % Saccade detection method (Default 'SVT')
params.saccadeThresh    = 1000;   % Saccade detection threshold (Default 1000)
params.saccadeLRPad     = 0.05;   % Padding in seconds to expand saccade candidates (Default 0.05)
params.minGoodChunk_len = 50;     % Minimum allowed gap between saccades (Default 50)

% File and folder paths (by default, all should be '')
params.exp_filepath    = '';  % Experiment file (.smr or .mat) path
params.cal_filepath    = '';  % Calibration file (.mat) path
params.save_folderpath = '';  % Target save folder path


%% BEGINNING OF SCRIPT SETUP

disp('--- RUNNING OKR NDD ANALYSIS ! ---');

% Run initial analysis setup GUI
[params, data] = run_okrNDD_setupGUI(params);

% If GUI was closed prematurely and not completed, throw an error to abort script
if ~all(params.setupAllCompleted)
    error('--- Analysis window was manually closed. ABORTING SCRIPT! ---');
else
    fprintf('\nAnalysis setup was successfully completed!\n')
end

% Extract the data from the GUI
nBlocks = data.nBlocks;
blockIds = data.blockIds;
blockTypes = data.blockTypes;
bssTimes = data.bssTimes;
cycleLength = data.cycleLength;
cycleTimes = data.cycleTimes;
chairvel_raw = data.chairvel_raw;
drumvel_raw = data.drumvel_raw;
hevel_raw = data.hevel_raw;
hepos_raw = data.hepos_raw;
saccadePad = data.saccadePad;
clear data;


%% BLOCK-LEVEL ANALYSIS

tic;
fprintf('\nStarting analysis:\n');
disp('    Analyzing data by blocks...');
tpbn = 0;  % Counter for keeping track of timepoint blocks
blocks = struct;
for ii = 1:nBlocks
    % Extract corresponding block segments for chair, drum, and eye data
    chairvel_raw_ii = chairvel_raw(blockIds(1,ii):blockIds(2,ii));
    drumvel_raw_ii  = drumvel_raw(blockIds(1,ii):blockIds(2,ii));
    hevel_raw_ii    = hevel_raw(blockIds(1,ii):blockIds(2,ii));
    hepos_raw_ii    = hepos_raw(blockIds(1,ii):blockIds(2,ii));
    % Note: In the real code:
    % If block is "VOR", the chair and eye traces are flipped
    % If block is "OKR", only the chair is flipped
    % I think this is to ensure that the stimulus always starts in the "positive" direction.
    % However, I feel that this obfuscates direction of "ipsi"/"contra"...

    % Define the length and time vectors for blocks
    blockLength = length(drumvel_raw_ii);
    blockTimes  = (0:blockLength-1) / params.fs;

    % Perform initial fit of the raw eye velocity data
    keep = abs(hevel_raw_ii) < 5*std(abs(hevel_raw_ii)) + mean(abs(hevel_raw_ii));
    [fit0,~] = calc_sineFit(params.exp_stimfreq, params.fs, keep, hevel_raw_ii);
    hevel_raw_fit = fit0.eyevel_fit;
    clear keep fit0;

    % Perform similar steps from original "desaccadeVel3" function:
    % 1. Remove transients from raw eye position segment (already performed earlier)
    % 2. Apply lowpass filter on raw eye position segment (hepos_ii)
    % 3. Applies SG filter on filtered eye position segment (hevel_ii)
    % 4. Calculate MSE of hevel_ii from hevel_raw_fit (calculated above)
    % 5. Defines (omitCenters) logical mask based on threshold value
    % 6. Defines (omitH) logical mask as omitCenters with padding (saccadeLRPad)
    hepos_ii = butterworthfilter(hepos_raw_ii, params.lowpassCutoff, params.fs, 'n',4);
    hevel_ii = movingslope(hepos_ii, params.filterWindow) * params.fs;

    % Apply saccade detection to the filtered eye velocity trace to get the necessary saccade masks
    [saccDilMask, saccMask, hevel_des_ii] = run_saccadeDetection(hevel_ii, ...
        hevel_raw_fit, ...
        params.saccadeThresh, ...
        saccadePad, ...
        params.minGoodChunk_len, ...
        params.saccadeMethod);

    % The fraction of NaN datapoints in the block corresponds to the fraction of the block that was masked as a saccade
    % This is calculated via the mean of the padded saccade mask omitH
    saccadeFrac = mean(saccDilMask);

    % Calculate sinusoidal fits of the desaccaded eye, chair, and drum velocity block data
    [bfit,stat] = calc_sineFit(params.exp_stimfreq, params.fs, ~saccDilMask, hevel_ii, chairvel_raw_ii, drumvel_raw_ii);

    % Compute the cycle fit for the desaccaded eye velocity trace
    eyevel_des_cyclefit = sin(2*pi*params.exp_stimfreq*cycleTimes + deg2rad(bfit.eyevel_rel_phase+180)) * bfit.eyevel_amp;

    % Calculate the start point of the first positive stimulus cycle
    startpt = max(1, round(mod(-bfit.ref_angle,360) / 360 * params.fs/params.exp_stimfreq));

    % Use segmentIntoCycles to get the estimated start indices of cycles
    cyclestart_ids_mat = segmentIntoCycles((1:blockLength), startpt, cycleLength);
    startpts = cyclestart_ids_mat(:,1);
    clear cyclestart_ids_mat;

    % Use segmentIntoCycles (previously VOR_breaktrace) to get the cycles of the data
    chairvel_cyclemat  = segmentIntoCycles(chairvel_raw_ii, startpt, cycleLength);
    drumvel_cyclemat   = segmentIntoCycles(drumvel_raw_ii, startpt, cycleLength);
    hevel_cyclemat     = segmentIntoCycles(hevel_ii, startpt, cycleLength); 
    hevel_des_cyclemat = segmentIntoCycles(hevel_des_ii, startpt, cycleLength); 
    saccmask_mat       = segmentIntoCycles(double(saccMask), startpt, cycleLength);
    mse_cyclemat       = segmentIntoCycles((hevel_ii-hevel_raw_fit).^2, startpt, cycleLength);

    % saccmask_mat is used to find the "bad" cycles that contain any NaN values
    badCycles = any(saccmask_mat, 2);

    % This is used to calculate how many "good" cycles there are
    nGoodCycles      = sum(~badCycles);
    [nTotalCycles,~] = size(hevel_cyclemat);

    % Calculate cycle-means and corresponding SEM estimates
    chairvel_cyclecvd = cycleVarianceDecomposition(chairvel_cyclemat);
    drumvel_cyclecvd = cycleVarianceDecomposition(drumvel_cyclemat);
    hevel_cyclecvd = cycleVarianceDecomposition(hevel_cyclemat);

    % A separate cycle mean is calculated using only "good" cycles 
    if nGoodCycles > 0
        hevel_good_cyclemat = hevel_cyclemat(~badCycles, :);
        hevel_good_cyclecvd = cycleVarianceDecomposition(hevel_good_cyclemat);
    else
        hevel_good_cyclemat = zeros(size(hevel_cyclemat));
        hevel_good_cyclecvd = cycleVarianceDecomposition(zeros(size(hevel_cyclemat(1,:))));
    end
    clear badCycles;

    % Calculate sinusoidal fit of the "good" eye, chair, and drum velocity cycle-means
    [cfit,~] = calc_sineFit(params.exp_stimfreq, params.fs, [], ...
        hevel_good_cyclecvd.cycleMean, chairvel_cyclecvd.cycleMean, drumvel_cyclecvd.cycleMean);

    % Populate results data for current block
    if ismember(ii, params.timepoint_ids)
        tpbn = tpbn + 1;
        blocks(ii).timePoint = sprintf('T%d', params.timepoint_groups(tpbn));
    else
        blocks(ii).timePoint = 'NA';
    end
    blocks(ii).blockType = char(blockTypes(ii));
    blocks(ii).blockNumber = ii;
    if bfit.chairvel_amp > 3
        blocks(ii).stimType = 'Chair';
        [good_rel_gain,good_rel_gain_SEM] = calc_eyeRelGainSEM(hevel_good_cyclecvd.amplitudeMean, ...
            hevel_good_cyclecvd.amplitudeSEM, ...
            chairvel_cyclecvd.amplitudeMean, ...
            chairvel_cyclecvd.amplitudeSEM);
        [good_rel_phase,good_rel_phase_SEM] = calc_eyeRelPhaseSEM(hevel_good_cyclecvd.phaseMean_deg, ...
            hevel_good_cyclecvd.phaseSEM_deg, ...
            chairvel_cyclecvd.phaseMean_deg, ...
            chairvel_cyclecvd.phaseSEM_deg);
    elseif bfit.drumvel_amp > 3
        blocks(ii).stimType = 'Drum';
        [good_rel_gain,good_rel_gain_SEM] = calc_eyeRelGainSEM(hevel_good_cyclecvd.amplitudeMean, ...
            hevel_good_cyclecvd.amplitudeSEM, ...
            drumvel_cyclecvd.amplitudeMean, ...
            drumvel_cyclecvd.amplitudeSEM);
        [good_rel_phase,good_rel_phase_SEM] = calc_eyeRelPhaseSEM(hevel_good_cyclecvd.phaseMean_deg, ...
            hevel_good_cyclecvd.phaseSEM_deg, ...
            drumvel_cyclecvd.phaseMean_deg, ...
            drumvel_cyclecvd.phaseSEM_deg);
    else
        blocks(ii).stimType = 'NoStim';
        good_rel_gain = nan;
        good_rel_gain_SEM = nan;
        good_rel_phase = nan;
        good_rel_phase_SEM = nan;
    end
    blocks(ii).stimFreq = params.exp_stimfreq;
    blocks(ii).hevel_amp = bfit.eyevel_amp;
    blocks(ii).hevel_phase = bfit.eyevel_phase;
    blocks(ii).hevel_rel_gain = bfit.eyevel_rel_gain;
    blocks(ii).hevel_rel_phase = bfit.eyevel_rel_phase;
    blocks(ii).good_amp = hevel_good_cyclecvd.amplitudeMean;
    blocks(ii).good_amp_SEM = hevel_good_cyclecvd.amplitudeSEM;
    blocks(ii).good_phase = hevel_good_cyclecvd.phaseMean_deg;
    blocks(ii).good_phase_SEM = hevel_good_cyclecvd.phaseSEM_deg;
    blocks(ii).good_rel_gain = good_rel_gain;
    blocks(ii).good_rel_gain_SEM = good_rel_gain_SEM;
    blocks(ii).good_rel_phase = good_rel_phase;
    blocks(ii).good_rel_phase_SEM = good_rel_phase_SEM;
    blocks(ii).saccadeFrac = saccadeFrac;
    blocks(ii).nGoodCycles = nGoodCycles;
    blocks(ii).nTotalCycles = nTotalCycles;

    % Block-specific data
    blocks(ii).startTime = bssTimes(1,ii);
    blocks(ii).endTime = bssTimes(2,ii);
    if bfit.chairvel_amp > 3
        blocks(ii).stimvel = chairvel_raw_ii;
    elseif bfit.drumvel_amp > 3
        blocks(ii).stimvel = drumvel_raw_ii;
    else
        blocks(ii).stimvel = zeros(size(blockTimes));
    end
    blocks(ii).hepos_raw = hepos_raw_ii;
    blocks(ii).hevel_raw = hevel_raw_ii;
    blocks(ii).hevel = hevel_ii;
    blocks(ii).hevel_des = hevel_des_ii;
    blocks(ii).hevel_raw_fit = hevel_raw_fit;
    blocks(ii).hevel_des_fit = bfit.eyevel_fit;

    % Cycle-specific data
    blocks(ii).cycleStartIds = startpts;
    blocks(ii).des_cyclefit = eyevel_des_cyclefit;
    blocks(ii).good_cyclefit = cfit.eyevel_rel_fit;
    if bfit.chairvel_amp > 3
        blocks(ii).stimvel_cyclemat = chairvel_cyclemat;
    elseif bfit.drumvel_amp > 3
        blocks(ii).stimvel_cyclemat = drumvel_cyclemat;
    else
        blocks(ii).stimvel_cyclemat = zeros(size(hevel_cyclemat));
    end
    blocks(ii).cyclemat = hevel_cyclemat;
    blocks(ii).des_cyclemat = hevel_des_cyclemat;
    blocks(ii).good_cyclemat = hevel_good_cyclemat;
    blocks(ii).mse_cyclemat = mse_cyclemat;

    % Additional metrics needed for Excel sheet
    blocks(ii).headamp = bfit.chairvel_amp;
    blocks(ii).headangle = bfit.chairvel_angle;
    blocks(ii).drumamp = bfit.drumvel_amp;
    blocks(ii).drumangle = bfit.drumvel_angle;
    blocks(ii).rsquare = stat(1);
    blocks(ii).pval = stat(3);
    blocks(ii).variance = sqrt(stat(4));

    % Include cycle variance decomposition results
    if bfit.chairvel_amp > 3
        blocks(ii).stimvel_cyclecvd = chairvel_cyclecvd;
        stimvel_amp = cfit.chairvel_amp;
        stimvel_phase = cfit.chairvel_angle;
    elseif bfit.drumvel_amp > 3
        blocks(ii).stimvel_cyclecvd = drumvel_cyclecvd;
        stimvel_amp = cfit.drumvel_amp;
        stimvel_phase = cfit.drumvel_angle;
    else
        blocks(ii).stimvel_cyclecvd = struct;
        stimvel_amp = nan;
        stimvel_phase = nan;
    end
    blocks(ii).hevel_cyclecvd = hevel_cyclecvd;
    blocks(ii).hevel_good_cyclecvd = hevel_good_cyclecvd;
    blocks(ii).stimvel_amp = blocks(ii).stimvel_cyclecvd.amplitudeMean;
    blocks(ii).stimvel_amp_SEM = blocks(ii).stimvel_cyclecvd.amplitudeSEM;
    blocks(ii).stimvel_phase = blocks(ii).stimvel_cyclecvd.phaseMean_deg;
    blocks(ii).stimvel_phase_SEM = blocks(ii).stimvel_cyclecvd.phaseSEM_deg;

    % Sanity check by ensuring fit results of sineFit and CVD agree
    check_valuesApproxEqual(cfit.eyevel_amp, blocks(ii).hevel_good_cyclecvd.amplitudeMean, '"Good" cycles amplitudes');
    check_valuesApproxEqual(cfit.eyevel_phase, blocks(ii).hevel_good_cyclecvd.phaseMean_deg, '"Good" cycles phases');
    check_valuesApproxEqual(stimvel_amp, blocks(ii).stimvel_cyclecvd.amplitudeMean, 'Stim cycles amplitudes');
    check_valuesApproxEqual(stimvel_phase, blocks(ii).stimvel_cyclecvd.phaseMean_deg, 'Stim cycles phases');
    check_valuesApproxEqual(cfit.eyevel_rel_gain, good_rel_gain, 'Relative gain');
    check_valuesApproxEqual(cfit.eyevel_rel_phase, good_rel_phase, 'Relative phase');
end

% Clean up workspace
clearvars('-except', 'analysis', 'blocks', 'cycleLength', 'cycleTimes', 'diffdata', 'params', 'timepoints');

% Generate Excel sheet of general block analysis results
if ~isequal(params.save_folderpath, 0) & isfield(params, 'save_filepath')
    disp('    Saving block data results to Excel sheet...');
    save_BlockAnalysisToXlsx(blocks, params);
end


%% TIMEPOINT GROUP ANALYSIS
disp('    Analyzing data grouped by timepoint...');

% Prepare necessary timepoint information
timepoints = struct;
tp_labels = string({blocks.timePoint});          % Timepoint labels of main rows
tp_labels = tp_labels(~strcmp(tp_labels,'NA'));  % Filter out all 'NA' rows
tp_unique = unique(tp_labels);                   % Unique timepoint labels
ntp = length(tp_unique);

for ii = 1:ntp
    tp_ii = tp_unique(ii);  % ii-th unique timepoint label
    tp_ii_ids = find(strcmp({blocks.timePoint}, tp_ii));
    tp_blocknumbers = strjoin(string([blocks(tp_ii_ids).blockNumber]), ', ');
    nids = length(tp_ii_ids);

    % Preallocate arrays for other cycle-means and cycle-medians
    timepoint_cyclelabels = [];  % Block number labels that corresponds to timepoint_cyclemat below
    good_cyclemean_mat   = zeros(nids, cycleLength);
    good_cyclemedian_mat = zeros(nids, cycleLength);

    % Loop through tp_ii_ids to collect the cycle-means and cycle-medians
    for jj = 1:nids
        id_jj = tp_ii_ids(jj);
        ngc_jj = size(blocks(id_jj).good_cyclemat, 1);
        timepoint_cyclelabels = [timepoint_cyclelabels; repmat(id_jj, ngc_jj, 1)];
        good_cyclemean_mat(jj,:)   = mean(blocks(id_jj).good_cyclemat, 1, 'omitnan');
        good_cyclemedian_mat(jj,:) = median(blocks(id_jj).good_cyclemat, 1, 'omitnan');
    end

    % Aggregate/pool all cycles from all blocks in the timepoint
    drumvel_cyclemat = vertcat(blocks(tp_ii_ids).stimvel_cyclemat);
    timepoint_cyclemat = vertcat(blocks(tp_ii_ids).good_cyclemat);

    % Perform Omnibus F-test on timepoint blocks to characterize block effects within the same timepoint
    res_ii = block_ftest(timepoint_cyclemat, timepoint_cyclelabels, 'Fs',params.fs, 'StimFreq',params.exp_stimfreq);

    % Calculate the mean and median of the pooled timepoint cycles
    nPooledCycles = size(timepoint_cyclemat, 1);
    pooled_drumvel_cyclemean_cvd = cycleVarianceDecomposition(drumvel_cyclemat);
    pooled_good_cyclemean_cvd = cycleVarianceDecomposition(timepoint_cyclemat);
    pooled_good_cyclemedian = median(timepoint_cyclemat, 1, 'omitnan');

    % Calculate mean of cycle-means and cycle-medians
    good_cyclemean_cvd = cycleVarianceDecomposition(good_cyclemean_mat);
    good_cyclemedian_cvd = cycleVarianceDecomposition(good_cyclemedian_mat);

    % Calculate nGoodCycles-weighted mean of "good" cycle-means
    nGoodCycles = [blocks(tp_ii_ids).nGoodCycles];
    good_nGoodCyclesWeighted_cyclemean = sum(good_cyclemean_mat.*nGoodCycles', 1) / sum(nGoodCycles');

    % Verify that the pooled good cycle-mean and nGoodCycles-weighted mean are equal
    check_valuesApproxEqual(pooled_good_cyclemean_cvd.cycleMean, ...
        good_nGoodCyclesWeighted_cyclemean, ...
        [tp_ii, ' cycle-means']);

    % Add to timepoint_results struct
    timepoints(ii).timePoint = char(tp_ii);
    timepoints(ii).blockType = blocks(tp_ii_ids(1)).blockType;
    timepoints(ii).blockNumbers = tp_blocknumbers;
    timepoints(ii).blockResultsIds = tp_ii_ids;

    % Calculate cycle metrics on each cycle average type
    pooled_good_cyclemean_cm = calc_cycleMetrics(pooled_good_cyclemean_cvd.cycleMean, pooled_drumvel_cyclemean_cvd.cycleMean);
    pooled_good_cyclemedian_cm = calc_cycleMetrics(pooled_good_cyclemedian, pooled_drumvel_cyclemean_cvd.cycleMean);
    good_cyclemean_cm = calc_cycleMetrics(good_cyclemean_cvd.cycleMean, pooled_drumvel_cyclemean_cvd.cycleMean);
    good_cyclemedian_cm = calc_cycleMetrics(good_cyclemedian_cvd.cycleMean, pooled_drumvel_cyclemean_cvd.cycleMean);
    nGoodCyclesWeighted_cyclemean_cm = calc_cycleMetrics(good_nGoodCyclesWeighted_cyclemean, pooled_drumvel_cyclemean_cvd.cycleMean);

    timepoints(ii).nPooledCycles = nPooledCycles;
    timepoints(ii).pooled_drumvel_cyclemean_cvd = pooled_drumvel_cyclemean_cvd;
    timepoints(ii).pooled_good_cyclemean_cvd = pooled_good_cyclemean_cvd;
    timepoints(ii).pooled_good_cyclemean_cm = pooled_good_cyclemean_cm;
    timepoints(ii).pooled_good_cyclemedian = pooled_good_cyclemedian;
    timepoints(ii).pooled_good_cyclemedian_cm = pooled_good_cyclemedian_cm;

    timepoints(ii).good_cyclemean_cvd = good_cyclemean_cvd;
    timepoints(ii).good_cyclemean_cm = good_cyclemean_cm;
    timepoints(ii).good_cyclemedian_cvd = good_cyclemedian_cvd;
    timepoints(ii).good_cyclemedian_cm = good_cyclemedian_cm;
    timepoints(ii).nGoodCyclesWeighted_cyclemean = good_nGoodCyclesWeighted_cyclemean;
    timepoints(ii).nGoodCyclesWeighted_cyclemean_cm = nGoodCyclesWeighted_cyclemean_cm;

    timepoints(ii).good_cyclemean_mat = good_cyclemean_mat;
    timepoints(ii).good_cyclemedian_mat = good_cyclemedian_mat;

    timepoints(ii).ftestFracSigPoints = res_ii.fracSigPoints;  % Fraction of cycle that were found significantly different
    timepoints(ii).ftestEtaSquared = res_ii.etaSquaredPartial;  % Eta Squared across the cycle
    timepoints(ii).ftestStats = res_ii.fStats;  % F-statistic across the cycle
end

% Clean up workspace
clearvars('-except', 'analysis', 'blocks', 'cycleLength', 'cycleTimes', 'diffdata', 'ntp', 'params', 'timepoints');


%% DIFFDATA ANALYSIS
disp('    Conducting diffdata analysis...');

% Conduct the diffdata analysis for permutations of timepoint pairs
diffdata = struct;
tp_pair_ids = get_diffdata_permutations(1:ntp);
npairs = size(tp_pair_ids, 1);
for ii = 1:npairs
    preidx = tp_pair_ids(ii,1);
    postidx = tp_pair_ids(ii,2);
    timePointDiff_label = sprintf('%s-%s', timepoints(postidx).timePoint, timepoints(preidx).timePoint);

    % Set label for ii-th timepoint pair
    diffdata(ii).timePointDiff = timePointDiff_label;
    diffdata(ii).timePointPrePostIds = [preidx, postidx];

    % Calculate post-pre for the various timepoint group cycle averages
    diffdata(ii).good_cyclemean_diff = timepoints(postidx).good_cyclemean_cvd.cycleMean - timepoints(preidx).good_cyclemean_cvd.cycleMean;
    diffdata(ii).good_cyclemean_diff_SEM = calc_cyclemeanDiffSEM(timepoints(preidx).good_cyclemean_cvd.cycleSEM, timepoints(postidx).good_cyclemean_cvd.cycleSEM);
    diffdata(ii).good_cyclemean_diff_cm = calc_cycleMetrics(diffdata(ii).good_cyclemean_diff, timepoints(preidx).pooled_drumvel_cyclemean_cvd.cycleMean);
    diffdata(ii).nGoodCyclesWeighted_cyclemean_diff = timepoints(postidx).nGoodCyclesWeighted_cyclemean - timepoints(preidx).nGoodCyclesWeighted_cyclemean;
    diffdata(ii).nGoodCyclesWeighted_cyclemean_diff_cm = calc_cycleMetrics(diffdata(ii).nGoodCyclesWeighted_cyclemean_diff, timepoints(preidx).pooled_drumvel_cyclemean_cvd.cycleMean);
    diffdata(ii).good_cyclemedianmean_diff = timepoints(postidx).good_cyclemedian_cvd.cycleMean - timepoints(preidx).good_cyclemedian_cvd.cycleMean;
    diffdata(ii).good_cyclemedianmean_diff_cm = calc_cycleMetrics(diffdata(ii).good_cyclemedianmean_diff, timepoints(preidx).pooled_drumvel_cyclemean_cvd.cycleMean);
    diffdata(ii).pooled_good_cyclemean_diff = timepoints(postidx).pooled_good_cyclemean_cvd.cycleMean - timepoints(preidx).pooled_good_cyclemean_cvd.cycleMean;
    diffdata(ii).pooled_good_cyclemean_diff_cm = calc_cycleMetrics(diffdata(ii).pooled_good_cyclemean_diff, timepoints(preidx).pooled_drumvel_cyclemean_cvd.cycleMean);
    diffdata(ii).pooled_good_cyclemean_diff_SEM = calc_cyclemeanDiffSEM(timepoints(preidx).pooled_good_cyclemean_cvd.cycleSEM, timepoints(postidx).pooled_good_cyclemean_cvd.cycleSEM);
    diffdata(ii).pooled_good_cyclemedian_diff = timepoints(postidx).pooled_good_cyclemedian - timepoints(preidx).pooled_good_cyclemedian;
    diffdata(ii).pooled_good_cyclemedian_diff_cm = calc_cycleMetrics(diffdata(ii).pooled_good_cyclemedian_diff, timepoints(preidx).pooled_drumvel_cyclemean_cvd.cycleMean);
end

% Add sheets to existing Excel for diffdata analysis results
if ~isequal(params.save_folderpath, 0) & isfield(params, 'save_filepath')
    disp('    Adding diffdata analysis results to Excel file...');
    save_DiffDataAnalysisToXlsx(timepoints, diffdata, params);
end

% Clean up workspace
clearvars('-except', 'analysis', 'blocks', 'diffdata', 'params', 'timepoints');

disp('    Compiling analysis results...');
analysis.params     = params;
analysis.blocks     = blocks;
analysis.timepoints = timepoints;
analysis.diffdata   = diffdata;
clear blocks timepoints diffdata;

fprintf('\nGenerating plots of analysis results:\n');
fignum = plot_standard_subplots(analysis, params.save_folderpath);
fignum = plot_okr_NDD_analysis(analysis, params.save_folderpath, fignum);

temp_analysis = analysis;
if ~isequal(params.save_folderpath, 0) & isfield(params, 'save_filepath')
    fprintf('    Saving analysis parameters to: %s\n', params.save_filepath); 
    clear analysis;
    analysis.params = temp_analysis.params;
    save(params.save_filepath, 'analysis');
    clear analysis;
    analysis = temp_analysis;
end

% Final cleanup
clearvars('-except', 'analysis');

fprintf('\n---ANALYSIS COMPLETED!---\n\n');
msgbox('Analysis complete!', 'Done');
toc


%% HELPER FUNCTIONS
function check_valuesApproxEqual(a, b, label, tol)
    if ~exist('tol', 'var')
        tol = 5e-7;
    end
    if any(isnan([a,b]))
        warning('NaN values for %s cannot be compared: a = %g vs b = %g', label, a, b);
    end
    status = all(abs(a-b) < tol);
    if ~status
        warning('Values for %s do not agree: a = %g vs b = %g', label, a, b);
    end
end

function diffSEM = calc_cyclemeanDiffSEM(cyclesem1, cyclesem2)
    cs1_sq = cyclesem1.^2;
    cs2_sq = cyclesem2.^2;
    diffSEM = sqrt(cs1_sq + cs2_sq);
end

function [gain,gainSEM] = calc_eyeRelGainSEM(eyeAmp, eyeSEM, stimAmp, stimSEM)
    gain = eyeAmp / stimAmp;
    a_sq = (eyeSEM / eyeAmp)^2;
    b_sq = (stimSEM / stimAmp)^2;
    gainSEM = gain * sqrt(a_sq + b_sq);
end

function [phase,phaseSEM] = calc_eyeRelPhaseSEM(eyePhase, eyeSEM, stimPhase, stimSEM)
    phase = (eyePhase - stimPhase);
    phase = mod(phase, 360) - 180;
    phaseSEM = sqrt((eyeSEM^2) + (stimSEM^2));
end