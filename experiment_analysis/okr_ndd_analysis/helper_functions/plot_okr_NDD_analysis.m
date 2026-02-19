function figure_number = plot_okr_NDD_analysis(analysis, savepath, figure_number)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
close all;

if ~exist('analysis', 'var')
    [filename,savepath] = uigetfile('*.mat', 'Select okr "NDD" analysis file', pwd);
    load(fullfile(savepath, filename), 'analysis');
    clear filename;
end

% Ensure savepath is a valid value (either 0 or a file path)
if ~exist('savepath', 'var')
    if ~isequal(analysis.params.save_folderpath,0) && isfolder(analysis.params.save_folderpath)
        savepath = uigetdir(analysis.params.save_folderpath, 'Select folder to save results in');
    else
        savepath = uigetdir(pwd, 'Select folder to save results in');
    end
end

if isequal(savepath , 0)
    warning('Save path was not provided. Figures will not be saved!');
end

if ~exist('figure_number','var') || isempty(figure_number) || ~isnumeric(figure_number)
    figure_number = 0;
end

% Load necessary variables and structs from analysis data
fs = analysis.params.fs;
params = analysis.params;
blocks = analysis.blocks;
timepoints = analysis.timepoints;
diffdata = analysis.diffdata;
clear analysis;

% Initialize length and time vectors for cycles (always the same)
cycleLength = round(params.fs / params.exp_stimfreq);
cycleTimes  = (1:cycleLength) / params.fs;

unique_timepoints = unique(string({blocks.timePoint}));
unique_timepoints = unique_timepoints(~strcmp(unique_timepoints,'NA'));  
ntimepoints = length(unique_timepoints);
plot_number = 0;
filename = {};


%% Figure: Blocks (Saccades and Cycle Segmentation)
for ii = 1:ntimepoints
    figure_number = figure_number + 1;
    plot_number = plot_number + 1;
    tp_ii = unique_timepoints(ii);
    tp_ids = timepoints(ii).blockResultsIds;
    br_ii = blocks(tp_ids);

    fig(plot_number) = figure('WindowState','maximized');
    h(plot_number) = tiledlayout(fig(plot_number), length(tp_ids), 1, 'TileSpacing','compact', 'Padding','compact');
    title(h(plot_number), sprintf("%s: %s BLOCKS (SACCADES AND CYCLE SEGMENTATION)", params.exp_filename, tp_ii), 'FontSize',14, 'Interpreter','none');
    filename{end+1} = sprintf('%02d_%s_blocks.fig', figure_number, tp_ii);

    ylimits = [-80, 80];
    nblocks = length(tp_ids);
    for jj = 1:nblocks
        btimes = (0:length(br_ii(jj).stimvel)-1) / fs;
        ax = nexttile;
        axesCustomToolbarButtons(ax);
        hold on;
        yline(0, 'Color',[0 0 0]+0.7, 'LineWidth',0.1, 'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
        plot(btimes, br_ii(jj).stimvel, '--k', 'DisplayName',sprintf('%s Velocity',br_ii(jj).stimType), 'HitTest','off', 'PickableParts','none'); 
        plot(btimes, br_ii(jj).hevel, '-r', 'DisplayName','Saccades', 'HitTest','off', 'PickableParts','none');
        plot(btimes, br_ii(jj).hevel_des, 'b', 'DisplayName','Filtered & Desaccaded Eye Velocity', 'HitTest','off', 'PickableParts','none'); 
        plot(btimes, br_ii(jj).hevel_des_fit, '-k', 'DisplayName','Fit of Eye Velocity', 'HitTest','off', 'PickableParts','none');
        % Preprocess line segments to improve performance of figures
        xpts = btimes(br_ii(jj).cycleStartIds);
        xdata = [xpts(:)'; xpts(:)'; nan(1,length(xpts))];
        ydata = repmat([ylimits(1)+0.25*abs(ylimits(1)); ylimits(2)-0.25*abs(ylimits(2)); NaN], 1, length(xpts));
        plot(xdata(:), ydata(:), 'Color','m', 'LineWidth',0.2, 'DisplayName','Computed Cycle Segments', 'HitTest','off', 'PickableParts','none');
        hold off; 
        title(sprintf('Block %d: %s (Gain = %.3f | Good Cycles = %d of %d | Saccade Fraction = %.3f)', ...
            br_ii(jj).blockNumber, br_ii(jj).blockType, br_ii(jj).hevel_rel_gain, br_ii(jj).nGoodCycles, br_ii(jj).nTotalCycles, br_ii(jj).saccadeFrac));
        xlim([btimes(1), btimes(end)]);
        ylim(ylimits);
        ylabel('Velocity (deg/s)');
        legend;
        box on;
    end

    if ~isequal(savepath, 0) && isfolder(savepath)
        savefig(fig(plot_number), fullfile(savepath, filename{end}));
        fprintf('    Figure %d saved as: %s\n', figure_number, filename{end});
    end
end


%% Figure: Cycles (Saccades and Cycle Segmentation)

for ii = 1:ntimepoints
    figure_number = figure_number + 1;
    plot_number = plot_number + 1;
    tp_ii = unique_timepoints(ii);
    tp_ids = timepoints(ii).blockResultsIds;
    br_ii = blocks(tp_ids);

    fig(plot_number) = figure('WindowState','maximized');
    h(plot_number) = tiledlayout(fig(plot_number), length(tp_ids), 3, 'TileSpacing','compact', 'Padding','compact');
    title(h(plot_number), sprintf("%s: %s CYCLES (GOOD CYCLE-MEANS VS CYCLE-MEDIANS | SACCADE THRESHOLD COMPARISONS)", params.exp_filename, tp_ii), 'FontSize',14, 'Interpreter','none');
    filename{end+1} = sprintf('%02d_%s_cycles.fig', figure_number, tp_ii);

    minval = 99;
    maxval = -99;
    pltnum = 0;
    nblocks = length(tp_ids);
    for jj = 1:nblocks
        [stimin, stimax] = bounds(br_ii(jj).stimvel_cyclecvd.cycleMean);
        if stimin < minval, minval = stimin; end
        if stimax > maxval, maxval = stimax; end

        pltnum = pltnum + 1;
        axs(pltnum) = nexttile;
        subplot_description = ['Compares the cycle-mean and cycle-median of "good" cycles. ', ...
            'Individual "good" cycles and stimulus velocity cycle-mean are overlayed in the background.'];
        axesCustomToolbarButtons(axs(pltnum),[],subplot_description);
        hold on;
        % Similarly preprocess lines to improve performance
        [nCycles, ~] = size(br_ii(jj).good_cyclemat);
        cycleData = [br_ii(jj).good_cyclemat'; nan(1, nCycles)];
        timeData = [repmat(cycleTimes(:), 1, nCycles); nan(1, nCycles)];
        plot(timeData(:), cycleData(:), 'Color',[1 0.7 1], 'LineWidth',0.1, 'DisplayName','Individual Eye Velocity Cycles', 'HitTest','off', 'PickableParts','none');
        yline(0, 'Color',[0 0 0]+0.7, 'LineWidth',0.1, 'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
        plot(cycleTimes, br_ii(jj).stimvel_cyclecvd.cycleMean, '--k', 'LineWidth',2, 'DisplayName',sprintf('%s Velocity Cycle-Mean',br_ii(jj).stimType), 'HitTest','off', 'PickableParts','none'); 
        plot(cycleTimes, br_ii(jj).good_cyclefit, '-k', 'LineWidth',2, 'DisplayName','Fit of Eye Velocity', 'HitTest','off', 'PickableParts','none');
        plot(cycleTimes, timepoints(ii).good_cyclemean_mat(jj,:), 'b', 'LineWidth',3, 'DisplayName','Eye Velocity Cycle-Mean', 'HitTest','off', 'PickableParts','none');
        plot(cycleTimes, timepoints(ii).good_cyclemedian_mat(jj,:), '--r', 'LineWidth',2, 'DisplayName','Eye Velocity Cycle-Median', 'HitTest','off', 'PickableParts','none');
        hold off;
        title(sprintf('Block %d: %s "Good" Cycle-Mean vs Cycle-Median (nGoodCycles = %d of %d | Gain = %.3f)', ...
            br_ii(jj).blockNumber, br_ii(jj).blockType, br_ii(jj).nGoodCycles, br_ii(jj).nTotalCycles, br_ii(jj).good_rel_gain));
        xlim([cycleTimes(1), cycleTimes(end)]);
        legend;
        box on;
        [cmin,cmax] = bounds(cycleData, 'all');
        if cmin < minval, minval = cmin; end
        if cmax > maxval, maxval = cmax; end

        % Build the shaded patch for SEM (upper then lower bound, traced as a polygon)
        upper = br_ii(jj).hevel_good_cyclecvd.cycleMean + br_ii(jj).hevel_good_cyclecvd.cycleSEM;
        lower = br_ii(jj).hevel_good_cyclecvd.cycleMean - br_ii(jj).hevel_good_cyclecvd.cycleSEM;
        semMean = mean(br_ii(jj).hevel_good_cyclecvd.cycleSEM);
        semStd = std(br_ii(jj).hevel_good_cyclecvd.cycleSEM);

        pltnum = pltnum + 1;
        axs(pltnum) = nexttile;
        subplot_description = ['Shows the eye velocity cycle-mean and corresponding SEM. ', ...
            'The mean and std of the SEM is also provided in the subplot title.'];
        axesCustomToolbarButtons(axs(pltnum),[],subplot_description);
        hold on;
        fill([cycleTimes, fliplr(cycleTimes)], [upper, fliplr(lower)], [1 0.7 1], 'FaceAlpha', 0.6, 'EdgeColor', 'none', ...
            'DisplayName','Cycle-Mean SEM', 'HitTest','off', 'PickableParts','none');
        yline(0, 'Color',[0 0 0]+0.7, 'LineWidth',0.1, 'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
        plot(cycleTimes, br_ii(jj).stimvel_cyclecvd.cycleMean, '--k', 'LineWidth',2, 'DisplayName',sprintf('%s Velocity Cycle-Mean',br_ii(jj).stimType), 'HitTest','off', 'PickableParts','none'); 
        plot(cycleTimes, br_ii(jj).hevel_good_cyclecvd.fitMean, '-k', 'LineWidth',2, 'DisplayName','Fit of Eye Velocity', 'HitTest','off', 'PickableParts','none');
        plot(cycleTimes, br_ii(jj).hevel_good_cyclecvd.cycleMean, '-b', 'LineWidth',2.5, 'DisplayName','Eye Velocity Cycle-Mean', 'HitTest','off', 'PickableParts','none');
        hold off;
        title(sprintf('Block %d: %s "Good" Cycle-Mean & SEM (mean = %.3f, std = %.3f)', ...
            br_ii(jj).blockNumber, br_ii(jj).blockType, semMean, semStd));
        xlim([cycleTimes(1), cycleTimes(end)]);
        legend;
        box on;

        ampVals = [br_ii(jj).stimvel_amp, br_ii(jj).hevel_cyclecvd.amplitudeMean, br_ii(jj).good_amp];
        ampSEMs = [br_ii(jj).stimvel_amp_SEM,  br_ii(jj).hevel_cyclecvd.amplitudeSEM, br_ii(jj).good_amp_SEM];
        phaseVals = [br_ii(jj).stimvel_phase,  br_ii(jj).hevel_cyclecvd.phaseMean_deg, br_ii(jj).good_phase];
        phaseSEMs = [br_ii(jj).stimvel_phase_SEM,  br_ii(jj).hevel_cyclecvd.phaseSEM_deg, br_ii(jj).good_phase_SEM];
        catLbls = {'Stimulus', 'AllCycles', 'GoodCyclesOnly'};
        ampLbls = arrayfun(@(m, s) sprintf('AmpMean = %.2f\n(AmpSEM = %.2f)', m, s), ...
                   ampVals, ampSEMs, 'UniformOutput', false);
        phsLbls = arrayfun(@(m, s) sprintf('PhaseMean = %.2f\n(PhaseSEM = %.2f)', m, s), ...
                   phaseVals, phaseSEMs, 'UniformOutput', false);

        maxAmpSEM = max(abs(ampSEMs));
        maxAmpSEM = maxAmpSEM + (0.5 * maxAmpSEM);
        maxPhsSEM = max(abs(phaseSEMs));
        maxPhsSEM = maxPhsSEM + (0.5 * maxPhsSEM);
        xVals = 1:length(catLbls);
        colororder({'b','r'});

        ax = nexttile;
        subplot_description = ['Compares the amplitude and phase estimates from the fits of the stimulus velocity, all eye velocity, and "good" eye velocity cycles.', ...
            'This can be used to simply assess the performance of the preprocessing on the accuracy of these estimates.'];
        axesCustomToolbarButtons(ax,[],subplot_description);
        hold on;
        yyaxis left;
        errorbar(xVals, ampVals, ampSEMs, 'o', 'Color','b');
        text(xVals+0.1, ampVals, ampLbls, 'Color','b', 'FontSize',10, 'HorizontalAlignment','left');
        ylim([min(ampVals)-maxAmpSEM, max(ampVals)+maxAmpSEM]);
        yyaxis right;
        errorbar(xVals, phaseVals, phaseSEMs, 'o', 'Color','r');
        text(xVals+0.1, phaseVals, phsLbls, 'Color','r', 'FontSize',10, 'HorizontalAlignment','left');
        ylim([min(phaseVals)-maxPhsSEM, max(phaseVals)+maxPhsSEM]);
        yyaxis left;
        hold off;
        title(sprintf('Block %d: %s Amplitude & Phase Estimate Comparisons', br_ii(jj).blockNumber, br_ii(jj).blockType));
        xlim([0.5, length(ampVals)+1.25]);
        xticks(xVals);         % Set tick positions
        xticklabels(catLbls);  % Set tick labels
        grid on;
        set(ax, 'GridAlpha',0.05);
        box on;
    end
    set(axs, 'YLim', [minval, maxval]);

    if ~isequal(savepath, 0) && isfolder(savepath)
        savefig(fig(plot_number), fullfile(savepath, filename{end}));
        fprintf('    Figure %d saved as: %s\n', figure_number, filename{end});
    end
end
clear ii jj kk br_ii calcCycleLimit cmin cmax cycleData cycleFracs cycleLimit cycleLimitThresh;
clear cycleSubMat_mean cycleSubMat_median linename minval maxval mseMaxValues nblocks nCycles nMinGoodCycles;
clear attempt rankIds rowNumber subplot_description timeData tp_ids tp_ii;


%% Figure: Mean of Cycle-Means Timepoint Differences
nTimePointDiffs = length(diffdata);
for ii = 1:nTimePointDiffs
    figure_number = figure_number + 1;
    plot_number = plot_number + 1;
    tpd_ii = diffdata(ii).timePointDiff;
    tr_ids = diffdata(ii).timePointPrePostIds;
    tr_ii = timepoints(tr_ids);
    cycleTimes = (1:length(tr_ii(1).pooled_drumvel_cyclemean_cvd.cycleMean)) / fs;
    halfCycle = round(length(cycleTimes) / 2);

    fig(plot_number) = figure('WindowState','maximized'); %'Units','normalized', 'Position',[0.01 0.2 0.95 0.6]
    h(plot_number) = tiledlayout(fig(plot_number), 7, 3, 'TileSpacing','compact', 'Padding','compact');
    title(h(plot_number), sprintf("%s: %s Pooled Cycle-Means and Cycle Difference", params.exp_filename, tpd_ii), 'FontSize',14, 'Interpreter','none');
    filename{end+1} = sprintf('%02d_%s_pooled_cyclediffs.fig', figure_number, tpd_ii);

    % Repeat first row plots for the pooled averages
    minval = 99;
    maxval = -99;
    npairs = length(tr_ii);
    for jj = 1:npairs
        [stimin, stimax] = bounds(tr_ii(jj).pooled_drumvel_cyclemean_cvd.cycleMean);
        if stimin < minval, minval = stimin; end
        if stimax > maxval, maxval = stimax; end
        if isequal(jj, 1)
            color1 = '-b';
            color2 = [0 0 0.4];
        else
            color1 = '-r';
            color2 = [0.4 0 0];
        end

        resjj = tr_ii(jj).pooled_good_cyclemean_cm;

        meanCen = resjj.eye.centroidMean(:);
        meanCenAdj = (meanCen + [0; halfCycle]) / fs;
        meanLag = resjj.centroidMeanDiff(:);
        meanSkew = resjj.eye.skewFromCentroidMean(:);

        medCen = resjj.eye.centroidMedian(:);
        medCenAdj = (medCen + [0; halfCycle]) / fs;
        medLag = resjj.centroidMedianDiff(:);
        medSkew = resjj.eye.skewFromCentroidMedian(:);

        bowleySkew = resjj.eye.skewQuantile(:);
        meanmedSkew = resjj.eye.skewMeanMedian(:);

        stimMeanMedSkew = mean(resjj.stim.skewMeanMedian(:), 'omitnan');

        meanTxt = cellstr(compose("MEAN (dotted):\n    Centroid = %.1f ms\n    Lag = %.1f ms\n    Skew = %.3f\n    MeanMedSkew = %.3f", meanCen(:), meanLag(:), meanSkew(:), meanmedSkew(:)));
        medTxt = cellstr(compose("MEDIAN (dashed):\n    Centroid = %.1f ms\n    Lag = %.1f ms\n    Skew = %.3f\n    BowleySkew = %.3f", medCen(:), medLag(:), medSkew(:), bowleySkew(:)));

        ax = nexttile([5 1]);
        subplot_description = [string(sprintf('Compares the mean and median of all POOLED "good" cycles of timepoint %s. ', tr_ii(jj).timePoint)) ...
        + "Various metrics to characterize 'peak' timing, lag time, and skew were computed on the pooled cycle-mean.", ...
        resjj.metadata.description];
        axesCustomToolbarButtons(ax,[],subplot_description);
        hold on;
        xline(meanCenAdj, 'k', meanTxt, 'LineStyle',':', 'LineWidth',0.001, 'LabelVerticalAlignment','middle', 'LabelOrientation','horizontal', 'Color',color2, 'HandleVisibility','off');
        xline(medCenAdj, 'k', medTxt, 'LineStyle','--', 'LineWidth',0.001, 'LabelVerticalAlignment','bottom', 'LabelOrientation','horizontal', 'Color',color2, 'HandleVisibility','off');
        xline(meanCenAdj, color1, 'LineStyle',':', 'LineWidth',1.5, 'HandleVisibility','off');
        xline(medCenAdj, color1, 'LineStyle','--', 'LineWidth',1, 'HandleVisibility','off');
        [nCycles, ~] = size(tr_ii(jj).good_cyclemean_mat);
        cycleData = [tr_ii(jj).good_cyclemean_mat'; nan(1, nCycles)];
        timeData = [repmat(cycleTimes(:), 1, nCycles); nan(1, nCycles)];
        plot(timeData(:), cycleData(:), 'Color',[1 0.6 1], 'LineWidth',1, 'DisplayName','Individual "Good" Eye Velocity Cycle-Means', 'HitTest','off', 'PickableParts','none');
        yline(0, 'Color',[0 0 0]+0.7, 'LineWidth',0.1, 'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
        plot(cycleTimes, tr_ii(jj).pooled_drumvel_cyclemean_cvd.cycleMean, '--k', 'LineWidth',2, 'DisplayName',sprintf('Drum Velocity MoCM (Mean-Median Skew = %.3f)', stimMeanMedSkew), 'HitTest','off', 'PickableParts','none');
        plot(cycleTimes, tr_ii(jj).pooled_good_cyclemean_cvd.cycleMean, color1, 'LineWidth',2, 'DisplayName',sprintf('Pooled Eye Velocity Cycle-mean'), 'HitTest','off', 'PickableParts','none');
        plot(cycleTimes, tr_ii(jj).pooled_good_cyclemedian, 'Color',[0 0.9 0.9], 'LineStyle',':', 'LineWidth',2.5, 'DisplayName','Pooled Eye Velocity Cycle-median', 'HitTest','off', 'PickableParts','none');
        hold off;
        title(sprintf('%s Blocks %s: Mean vs Median of POOLED "Good" Cycles (F-test FSP = %.3f)', tr_ii(jj).timePoint, tr_ii(jj).blockNumbers, tr_ii(jj).ftestFracSigPoints));
        xlim([0, cycleTimes(end)]);
        xlabel('Time (s)');
        if isequal(jj,1), ylabel('Velocity (deg/s)'); end
        legend;
        box on;
        [cmin,cmax] = bounds(cycleData, 'all');
        if cmin < minval, minval = cmin; end
        if cmax > maxval, maxval = cmax; end
    end

    cyclemeanDiff = diffdata(ii).pooled_good_cyclemean_diff;
    resjj = diffdata(ii).pooled_good_cyclemean_diff_cm;

    meanCen = resjj.eye.centroidMean(:);
    meanCenAdj = (meanCen + [0; halfCycle]) / fs;
    meanLag = resjj.centroidMeanDiff(:);
    meanSkew = resjj.eye.skewFromCentroidMean(:);

    medCen = resjj.eye.centroidMedian(:);
    medCenAdj = (medCen + [0; halfCycle]) / fs;
    medLag = resjj.centroidMedianDiff(:);
    medSkew = resjj.eye.skewFromCentroidMedian(:);

    bowleySkew = resjj.eye.skewQuantile(:);
    meanmedSkew = resjj.eye.skewMeanMedian(:);

    meanTxt = cellstr(compose("MEAN (dotted):\n    Centroid = %.1f ms\n    Lag = %.1f ms\n    Skew = %.3f\n    MeanMedSkew = %.3f", meanCen(:), meanLag(:), meanSkew(:), meanmedSkew(:)));
    medTxt = cellstr(compose("MEDIAN (dashed):\n    Centroid = %.1f ms\n    Lag = %.1f ms\n    Skew = %.3f\n    BowleySkew = %.3f", medCen(:), medLag(:), medSkew(:), bowleySkew(:)));

    % Find the two tallest peaks in the cycle-mean difference data
    pk1 = [resjj.eye.peak1Val(1), resjj.eye.peak2Val(1)];
    locs1 = (round([resjj.eye.peak1TimeMs(1), resjj.eye.peak2TimeMs(1)]) + 1) / fs;
    pk2 = [resjj.eye.peak1Val(2), resjj.eye.peak2Val(2)];
    locs2 = (round([resjj.eye.peak1TimeMs(2), resjj.eye.peak2TimeMs(2)]) + 1 + round(halfCycle)) / fs;

    ax = nexttile([5 1]);
    subplot_description = [string(sprintf('Compares the %s (blue), %s (red), and the diff of means %s-%s (green) of pooled "good" cycle-means. ', tr_ii(1).timePoint, tr_ii(2).timePoint, tr_ii(2).timePoint, tr_ii(1).timePoint)) ...
        + string(sprintf('Various metrics to characterize "peak" timing, lag time, and skew were computed on the diff of means %s-%s:', tr_ii(2).timePoint, tr_ii(1).timePoint)), ...
        resjj.metadata.description];
    axesCustomToolbarButtons(ax,[],subplot_description);
    hold on;
    yline(0, 'Color',[0 0 0]+0.7, 'LineWidth',0.1, 'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
    xline(meanCenAdj, 'k', meanTxt, 'LineStyle',':', 'LineWidth',0.001, 'LabelVerticalAlignment','middle', 'LabelOrientation','horizontal', 'Color',[0 0.4 0], 'HandleVisibility','off');
    xline(medCenAdj, 'k', medTxt, 'LineStyle','--', 'LineWidth',0.001, 'LabelVerticalAlignment','bottom', 'LabelOrientation','horizontal', 'Color',[0 0.4 0], 'HandleVisibility','off');
    xline(meanCenAdj, 'g', 'LineStyle',':', 'LineWidth',1.5, 'HandleVisibility','off');
    xline(medCenAdj, 'g', 'LineStyle','--', 'LineWidth',1, 'HandleVisibility','off');
    plot(cycleTimes, tr_ii(1).pooled_drumvel_cyclemean_cvd.cycleMean, '--k', 'LineWidth',2, 'DisplayName','Drum Velocity MoCM', 'HitTest','off', 'PickableParts','none');
    plot(cycleTimes, tr_ii(1).pooled_good_cyclemean_cvd.cycleMean, '-b', 'LineWidth',2, 'DisplayName',sprintf('%s Pooled Cycle-mean',tr_ii(1).timePoint), 'HitTest','off', 'PickableParts','none');
    plot(cycleTimes, tr_ii(2).pooled_good_cyclemean_cvd.cycleMean, '-r', 'LineWidth',2, 'DisplayName',sprintf('%s Pooled Cycle-mean',tr_ii(2).timePoint), 'HitTest','off', 'PickableParts','none');
    plot(cycleTimes, cyclemeanDiff, '-g', 'LineWidth',2, 'DisplayName',sprintf('%s-%s Pooled Cycle-mean Diff',tr_ii(2).timePoint,tr_ii(1).timePoint), 'HitTest','off', 'PickableParts','none');  
    plot(locs1, pk1, 'm.', 'MarkerSize',20, 'DisplayName','Peaks');
    plot(locs2, pk2, 'm.', 'MarkerSize',20, 'HandleVisibility','off');
    text(locs1, pk1, cellstr(" "+string(round(locs1*fs)))+"ms", 'Color','m', 'VerticalAlignment','bottom');
    text(locs2, pk2, cellstr(" "+string(round(locs2*fs)))+"ms", 'Color','m', 'VerticalAlignment','top');
    hold off;
    title(sprintf('Timepoints %s vs %s vs %s-%s ("Diff of Means"): Pooled "Good" Cycle-mean', tr_ii(1).timePoint, tr_ii(2).timePoint, tr_ii(2).timePoint, tr_ii(1).timePoint));
    xlim([0, cycleTimes(end)]);
    legend;
    box on;

    ylims = [minval-0.2*abs(minval), maxval+0.2*abs(maxval)];
    set(findobj(h(plot_number),'Type','axes'), 'YLim', ylims);

    % Plot row corresponding to F-Test statistics of each block group
    for jj = 1:npairs
        % Plot row corresponding to F-Test statistics of each block group
        ax = nexttile([2 1]);
        subplot_description = ['Shows the Partial Eta-Squared effect size across the cycle. ', ...
            'This measures the proportion of total variance at each time point that is explained by block differences. ', ...
            'Values range from 0 to 1, where higher values indicate a stronger effect independent of sample size.'];
        axesCustomToolbarButtons(ax,[],subplot_description);
        hold on;
        area(cycleTimes, tr_ii(jj).ftestEtaSquared, 'HitTest','off', 'PickableParts','none');
        hold off;
        xlim([0, cycleTimes(end)]);
        ylim([0 1]);
        xlabel('Time (s)');
        title('Proportion of Total Variance');
        grid on;
        set(ax, 'GridAlpha',0.05);
        box on;
    end

    % Empty extra axis since difference doesn't have an F-test
    ax = nexttile([1 1]);
    ax.Visible = 'off';

    if ~isequal(savepath, 0) && isfolder(savepath)
        savefig(fig(plot_number), fullfile(savepath, filename{end}));
        fprintf('    Figure %d saved as: %s\n', figure_number, filename{end});
    end
end
clear tp_ii tp_ids br_ii minval maxval nblocks stimin stimax subplot_description;


%% Figure: Cycle Characterization of All Blocks  (TODO)
figure_number = figure_number + 1;
plot_number = plot_number + 1;

fig(plot_number) = figure('WindowState','maximized');
h(plot_number) = tiledlayout(fig(plot_number), 2, 3, 'TileSpacing','compact', 'Padding','compact');
title(h(plot_number), sprintf('%s: All Block "Good" Cycle-Means Analysis', params.exp_filename), 'FontSize',14, 'Interpreter','none');
filename{end+1} = sprintf('%02d_all_block_goodcycles_analysis.fig', figure_number);

nblocks = length(blocks);
ncycles = length(blocks(1).stimvel_cyclecvd.cycleMean);
good_cyclemeans = zeros(nblocks-2, ncycles);
stimvel_cyclemeans = zeros(nblocks-2, ncycles);
for ii = 1:nblocks-2
    good_cyclemeans(ii,:) = mean(blocks(ii+1).good_cyclemat, 1, 'omitnan');
    stimvel_cyclemeans(ii,:) = blocks(ii+1).stimvel_cyclecvd.cycleMean;
end

res = calc_cycleMetrics(good_cyclemeans, stimvel_cyclemeans);
blocknumbers = 2:nblocks-1;

ax = nexttile;
subplot_description = ['Measure of the temporal shift between response and stimulus centers of mass.', ...
    ' Positive values indicate the eye response lags the stimulus; negative values indicate the eye response leads.'];
axesCustomToolbarButtons(ax,[],subplot_description);
hold on;
yline(0, 'Color','k', 'LineStyle','--', 'LineWidth',1, 'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
scatter(blocknumbers, res.centroidMeanDiff(:,1), 25, 'b', 'filled', 'DisplayName','NT Half-Cycles');
scatter(blocknumbers, res.centroidMeanDiff(:,2), 25, 'r', 'DisplayName','TN Half-Cycles');
hold off;
xlim(ax, [1, blocknumbers(end)+1]);
title(ax, 'Centroid Mean Difference (Lead/Lag)', 'FontSize',14);
legend(ax);
grid on;

ax = nexttile;
subplot_description = ['Measures where the center-of-mass falls relative to the temporal midpoint of the half-cycle.', ...
    ' Positive or negative values correspond to early or late mass concentration.'];
axesCustomToolbarButtons(ax,[],subplot_description);
hold on;
yline(0, 'Color','k', 'LineStyle','--', 'LineWidth',1, 'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
scatter(blocknumbers, res.eye.skewFromCentroidMean(:,1), 25, 'b', 'filled');
scatter(blocknumbers, res.eye.skewFromCentroidMean(:,2), 25, 'r');
hold off;
xlim(ax, [1, blocknumbers(end)+1]);
title(ax, 'Skew via Centroid Mean', 'FontSize',14);
grid on;

ax = nexttile;
subplot_description = ['Quantifies how much tails pull the mean away from the median, normalized by the spread of the distribution.', ...
    ' Positive = mean is later than median, Negative = mean is earlier than median.'];
axesCustomToolbarButtons(ax,[],subplot_description);
hold on;
yline(0, 'Color','k', 'LineStyle','--', 'LineWidth',1, 'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
scatter(blocknumbers, res.eye.skewMeanMedian(:,1), 25, 'b', 'filled');
scatter(blocknumbers, res.eye.skewMeanMedian(:,2), 25, 'r');
hold off;
xlim(ax, [1, blocknumbers(end)+1]);
title(ax, 'Mean-Median Skew', 'FontSize',14);
grid on;

ax = nexttile;
subplot_description = ['Same as for centroid mean but using the median-based centroid. More robust to waveform distortions.', ...
    ' Positive = lag, Negative = lead.'];
axesCustomToolbarButtons(ax,[],subplot_description);
hold on;
yline(0, 'Color','k', 'LineStyle','--', 'LineWidth',1, 'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
scatter(blocknumbers, res.centroidMedianDiff(:,1), 25, 'b', 'filled');
scatter(blocknumbers, res.centroidMedianDiff(:,2), 25, 'r');
hold off;
xlim(ax, [1, blocknumbers(end)+1]);
xlabel(ax, 'Block Number');
title(ax, 'Centroid Median Difference (Lead/Lag)', 'FontSize',14);
grid on;

ax = nexttile;
subplot_description = ['Same as for centroid mean but using the median-based centroid.', ...
    ' Positive = median early, Negative = median late.'];
axesCustomToolbarButtons(ax,[],subplot_description);
hold on;
yline(0, 'Color','k', 'LineStyle','--', 'LineWidth',1, 'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
scatter(blocknumbers, res.eye.skewFromCentroidMedian(:,1), 25, 'b', 'filled');
scatter(blocknumbers, res.eye.skewFromCentroidMedian(:,2), 25, 'r');
hold off;
xlim(ax, [1, blocknumbers(end)+1]);
xlabel(ax, 'Block Number');
title(ax, 'Skew via Centroid Median', 'FontSize',14);
grid on;

ax = nexttile;
subplot_description = ['Measures asymmetry in the spread around the median.', ...
    ' Positive = late/right tail, Negative = early/left tail.'];
axesCustomToolbarButtons(ax,[],subplot_description);
hold on;
yline(0, 'Color','k', 'LineStyle','--', 'LineWidth',1, 'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
scatter(blocknumbers, res.eye.skewQuantile(:,1), 25, 'b', 'filled');
scatter(blocknumbers, res.eye.skewQuantile(:,2), 25, 'r');
hold off;
xlim(ax, [1, blocknumbers(end)+1]);
xlabel(ax, 'Block Number');
title(ax, 'Bowley Skewness (Quantile Skew)', 'FontSize',14);
grid on;

if ~isequal(savepath, 0) && isfolder(savepath)
    savefig(fig(plot_number), fullfile(savepath, filename{end}));
    fprintf('    Figure %d saved as: %s\n', figure_number, filename{end});
end


% Re-arrange the opened figures
for ii = length(fig):-1:1
    figure(fig(ii));
end


end