function figure_number = plot_standard_subplots(analysis, savepath, figure_number)
%UNTITLED5 Summary of this function goes here
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

if ~exist('figure_number','var') || isempty(figure_number) || ~isnumeric(figure_number)
    figure_number = 0;
end

if isequal(savepath , 0)
    warning('Save path was not provided. Figures will not be saved!');
    return
end

% Load necessary variables and structs from analysis data
fs = analysis.params.fs;
stimFreq = analysis.params.exp_stimfreq;
params = analysis.params;
blocks = analysis.blocks;

% Preallocate arrays for other cycle-means and cycle-medians
nblocks = length(blocks);
cycleLength = round(fs / params.exp_stimfreq);
ctimes = (0:cycleLength-1) / fs;
clear analysis;

% Initialize figure for subplots1
figure_number = figure_number + 1;
fig1 = figure('Units','pixels', 'Position',[200 200 1200 1000], 'Visible','off'); 
filename1 = sprintf('%02d_subplots1.pdf', figure_number);
s1a_ylim = [-80, 80];

% Initialize figure for subplots2
figure_number = figure_number + 1;
fig2 = figure('Units','pixels', 'Position',[200 200 1200 1000], 'Visible','off');
filename2 = sprintf('%02d_subplots2.pdf', figure_number);

% Plot subplot rows for each experiment block
N = 4;  % Number of subplot rows per page
pageCount = 0;  % Page counter
warning('off', 'MATLAB:print:ContentTypeImageSuggested');
for ii = 1:nblocks
    % At the start of each page, create fresh layouts
    rowOnPage = mod(ii-1, N) + 1;  % 1-based row index within current page
    if rowOnPage == 1
        clf(fig1);
        h1 = tiledlayout(fig1, N, 9, 'Padding','compact', 'TileSpacing','compact');
        clf(fig2);
        h2 = tiledlayout(fig2, N, 3, 'Padding','compact', 'TileSpacing','compact');
    end
    
    % Get relevant data for this block
    b_ii = blocks(ii);
    btimes = (0:length(b_ii.hevel)-1) / fs;

    % Set ylimits based on stimulus velocity
    ylimits = double([round(min(b_ii.stimvel_cyclecvd.cycleMean)/10)*12, ...
        round(max(b_ii.stimvel_cyclecvd.cycleMean)/10)*12]);
    
    % Check for bad or missing stim
    if abs(ylimits(1)) <= 2
        ylimits(1) = -2;
        ylimits(2) = 2;
    end

    %% FIGURE 1 SUBPLOTS

    % --- Subplot-1A: Block, Fit, and Saccades ---------------------------------- %
    ax = nexttile(h1, [1,6]);
    hold(ax,'on');
    plot(ax, btimes, b_ii.hevel, 'k', 'LineWidth', .5, 'DisplayName','Saccades');
    plot(ax, btimes, b_ii.hevel_des, 'b', 'LineWidth', .5, 'DisplayName','Desaccade Eye Velocity');
    plot(ax, btimes, b_ii.hevel_des_fit,'r', 'LineWidth', .3, 'DisplayName','Sinusoidal Fit');
    hold(ax,'off');

    % Cosmetics
    xlim(ax, [0, btimes(end)]);
    ylim(ax, s1a_ylim);
    xticks(ax, linspace(0, max(xlim), 5));
    yticks(ax, linspace(s1a_ylim(1), s1a_ylim(2), 5));
    ylabel(ax, 'Velocity (deg/s)');
    title(ax, sprintf('%d min %s %gHz', round(b_ii.blockNumber-1), b_ii.blockType, b_ii.stimFreq));

    % Text displaying absolute time of segment start
    text(ax, 0, max(ylim)*1.15, ['@ ' num2str(round(b_ii.startTime, 2)), 's'], 'FontSize', 7);

    if isequal(rowOnPage,N) || isequal(ii,nblocks)
        xlabel(ax, 'Time (s)');
    end

    % --- Subplot-1B: Cycle and Fit ---------------------------------- %
    ax = nexttile(h1, [1,2]);
    hold(ax,'on');
    plot(ax, ctimes, b_ii.hevel_good_cyclecvd.cycleMean, 'b', 'DisplayName','"Good" Cycles-Only');
    plot(ax, ctimes, mean(b_ii.des_cyclemat, 1, 'omitnan'), 'g', 'DisplayName','All Desaccaded Cycles');
    plot(ax, ctimes, b_ii.des_cyclefit, 'r', 'DisplayName','Fit of All Desaccaded Cycles');
    plot(ax, ctimes, b_ii.stimvel_cyclecvd.cycleMean, 'k', 'DisplayName','Stimulus Velocity');
    yline(ax, 0, ':k', 'HandleVisibility','off');
    hold(ax,'off');

    % Cosmetics
    box(ax, 'off');
    ylim(ax, ylimits);
    xlim(ax, [0 round(ctimes(end))]);
    xticks(ax, linspace(0, max(xlim), 5));
    yticks(ax, linspace(min(ylim), max(ylim), 5));

    if isequal(rowOnPage,N) || isequal(ii,nblocks)
        xlabel(ax, 'Time (s)');
    end

    % --- Subplot-1C: Reference Information ---------------------------------- %
    ax = nexttile(h1, [1,1]);
    xtxt = min(xlim(ax)) + (0.1 * max(ylim(ax)));
    ytxt = max(ylim(ax));
    sep = (max(ylim(ax)) - 0.5*min(ylim(ax))) / 8;
    text(ax, xtxt, ytxt-1*sep, sprintf('Good Cycles: %d / %d', b_ii.nGoodCycles, b_ii.nTotalCycles), 'FontSize',7);
    text(ax, xtxt, ytxt-2*sep, sprintf('Rel Gain: %g', b_ii.hevel_rel_gain), 'FontSize',7);
    text(ax, xtxt, ytxt-3*sep, sprintf('Eye Amp: %.3f', b_ii.hevel_amp), 'FontSize',7);
    text(ax, xtxt, ytxt-4*sep, sprintf('Rel Phase: %.3f', b_ii.hevel_rel_phase), 'FontSize',7);
    text(ax, xtxt, ytxt-5*sep, sprintf('Stim: %s', b_ii.stimType), 'FontSize',7);
    text(ax, xtxt, ytxt-6*sep, sprintf('r^2: %g', b_ii.rsquare), 'FontSize',7, 'Interpreter','none');
    text(ax, xtxt, ytxt-7*sep, sprintf('sacFrac: %g', b_ii.saccadeFrac), 'FontSize',7);
    box(ax, 'off');
    grid(ax, 'off');
    axis(ax, 'off');

    %% FIGURE 2 SUBPLOTS

    % Skip plotting of blocks with 0 "good" cycles
    if (b_ii.nGoodCycles <= 0)
        ax = nexttile(h2); ax = nexttile(h2); ax = nexttile(h2);
    else
        % Preprocess lines to improve performance
        gcylims = quantile(blocks(ii).good_cyclemat, [0.001, 0.999], 'all');
        [nCycles,~] = size(b_ii.good_cyclemat);
        cycleData = [b_ii.good_cyclemat'; nan(1, nCycles)];
        timeData = [repmat(ctimes', 1, nCycles); nan(1, nCycles)];

        ax = nexttile(h2);
        if isequal(rowOnPage, 1), title(ax, 'Good Cycles and Cycle-mean'); end
        hold(ax, 'on');
        plot(ax, ctimes, b_ii.stimvel_cyclecvd.cycleMean, '--k', 'LineWidth', 1);
        plot(ax, timeData, cycleData, 'b', 'LineWidth', .1); 
        plot(ax, ctimes, b_ii.hevel_good_cyclecvd.cycleMean, 'k', 'LineWidth', 2);
        yline(ax, 0, ':k');
        hold(ax, 'off');

        % Cosmetics
        xlim(ax, [0 ctimes(end)]);
        ylim(ax, gcylims);
        ylabel(ax, sprintf('%d / %d', b_ii.nGoodCycles, b_ii.nTotalCycles), 'FontSize',12);
        if isequal(rowOnPage,N) || isequal(ii,nblocks)
            xlabel(ax, 'Time (s)');
        end
        box(ax, 'on');
        
        % Figure (B) Raster Plot of Residuals (Cycle - CycleFit)
        [nCycles2,~] = size(b_ii.cyclemat);
        mat_Residuals = (b_ii.cyclemat' - b_ii.good_cyclefit(:)) ./ b_ii.good_amp;
        badCycles = any(isnan(b_ii.des_cyclemat'), 1);
        mat_Residuals(:,badCycles) = NaN;
        alphaMap = double(~isnan(mat_Residuals'));
        
        ax = nexttile(h2);
        if isequal(rowOnPage, 1)
            title(ax, {'Good Cycle Residuals (Relative to Fit)', sprintf('%d min %s %gHz', round(b_ii.blockNumber-1), b_ii.blockType, b_ii.stimFreq)});
        else
            title(ax, sprintf('%d min %s %gHz', round(b_ii.blockNumber-1), b_ii.blockType, b_ii.stimFreq));
        end
        hold(ax, 'on');
        imagesc(ax, ctimes, 1:nCycles2, mat_Residuals', 'AlphaData',alphaMap);
        colormap(ax, 'jet');
        clim(ax, [min(mat_Residuals(:)), max(mat_Residuals(:))]);
        yyaxis right;
        plot(ax, ctimes, b_ii.good_cyclefit, '--k', 'LineWidth', 1);
        plot(ax, ctimes, b_ii.hevel_good_cyclecvd.cycleMean, '-k', 'LineWidth', 2);
        yline(ax, 0, ':k');
        ylim(ax, gcylims);
        yyaxis left;
        hold(ax, 'off');

        % Cosmetics
        xlim(ax, [0 ctimes(end)]);
        ylim(ax, [0 nCycles2+1]);
        if isequal(rowOnPage,N) || isequal(ii,nblocks)
            xlabel(ax, 'Time (s)');
        end
        ylabel(ax, 'Cycle Number', 'FontSize',12);
        box(ax, 'on');
        ax.YDir = 'reverse';

        % Figure (C) Power Spectrum of Filtered Position Traces
        hepos_raw_ii = b_ii.hepos_raw;
        hepos_raw_ii = hepos_raw_ii - mean(hepos_raw_ii);  % Remove DC offset
        [psd_mt, f_mt] = pmtm(hepos_raw_ii, 4, [], fs, 'DropLastTaper',true); 
        % --- 1/f slope estimate (fit log-log PSD between 0.1 and 10 Hz) -----
        idx_1f          = f_mt >= 0.1 & f_mt <= 10;
        % Exclude the 1 Hz VOR peak region for the fit
        idx_1f_nopeak   = idx_1f & ~(f_mt >= 0.8 & f_mt <= 1.2);
        log_f           = log10(f_mt(idx_1f_nopeak));
        log_psd         = log10(psd_mt(idx_1f_nopeak));
        p_1f            = polyfit(log_f, log_psd, 1);  % slope, intercept
        slope_1f        = p_1f(1);   % should be ~ -1 for pure 1/f noise, ~ -2 for 1/f^2

        ax = nexttile(h2);
        if isequal(rowOnPage, 1), title(ax, 'Multitaper Spectrum '); end
        hold(ax, 'on');
        plot(ax, f_mt, pow2db(psd_mt), 'b', 'LineWidth',1.2, 'DisplayName','Raw Eye Position', 'HitTest','off', 'PickableParts','none');
        xline(stimFreq, '--r', 'LineWidth',0.8, 'DisplayName','Stimulus frequency', 'HitTest','off', 'PickableParts','none');
        f_fit = logspace(log10(0.1), log10(params.lowpassCutoff), 100);
        plot(f_fit, 10*(polyval(p_1f, log10(f_fit))), 'k--', 'LineWidth',1, ...
            'DisplayName',sprintf('Log-log fit slope = %.2f',slope_1f), 'HitTest','off', 'PickableParts','none');
        hold(ax, 'off');

        % Cosmetics
        xlim(ax, [0.01 params.lowpassCutoff+(0.2*params.lowpassCutoff)]);
        if isequal(rowOnPage,N) || isequal(ii,nblocks)
            xlabel(ax, 'Frequency (Hz)');
        end
        legend(ax);
        box(ax, 'on');
    end

    %% EXPORT when page is full
    if rowOnPage == N
        pageCount = pageCount + 1;
        a1 = annotation(fig1, 'rectangle', [0 0 1 1], 'Color', 'w');
        a2 = annotation(fig2, 'rectangle', [0 0 1 1], 'Color', 'w');
        if pageCount == 1
            exportgraphics(fig1, fullfile(savepath, filename1), 'ContentType','vector');
            exportgraphics(fig2, fullfile(savepath, filename2), 'ContentType','vector');
        else
            exportgraphics(fig1, fullfile(savepath, filename1), 'ContentType','vector', 'Append',true);
            exportgraphics(fig2, fullfile(savepath, filename2), 'ContentType','vector', 'Append',true);
        end
        fprintf('    Page %d exported (blocks %d–%d).\n', pageCount, ii-rowOnPage+1, ii);
        delete(a1);
        delete(a2);
    end

end

% Handle final partial page: pad with empty hidden axes, then export
if mod(nblocks, N) ~= 0
    remainingRows = N - mod(nblocks, N);
    for jj = 1:remainingRows
        % Fig1: 9 tiles per row (6+2+1)
        ax = nexttile(h1, [1,6]); axis(ax,'off');
        ax = nexttile(h1, [1,2]); axis(ax,'off');
        ax = nexttile(h1, [1,1]); axis(ax,'off');
        % Fig2: 3 tiles per row
        ax = nexttile(h2); axis(ax,'off');
        ax = nexttile(h2); axis(ax,'off');
        ax = nexttile(h2); axis(ax,'off');
    end
    pageCount = pageCount + 1;
    a1 = annotation(fig1, 'rectangle', [0 0 1 1], 'Color', 'w');
    a2 = annotation(fig2, 'rectangle', [0 0 1 1], 'Color', 'w');
    if pageCount == 1
        exportgraphics(fig1, fullfile(savepath, filename1), 'ContentType','vector');
        exportgraphics(fig2, fullfile(savepath, filename2), 'ContentType','vector');
    else
        exportgraphics(fig1, fullfile(savepath, filename1), 'ContentType','vector', 'Append',true);
        exportgraphics(fig2, fullfile(savepath, filename2), 'ContentType','vector', 'Append',true);
    end
    delete(a1);
    delete(a2);
    fprintf('    Page %d exported (blocks %d–%d, padded).\n', pageCount, nblocks-mod(nblocks,N)+1, nblocks);
end
warning('on', 'MATLAB:print:ContentTypeImageSuggested');

end