function figure_number = plot_lpcgcComparison(analysis, T, D, savepath, figure_number)
% PLOT_LPCGCCOMPARISON  Draws the 5-tab lowpass-cutoff × saccade-threshold
% comparison figure (MAIN / FILTER / THRESH / ROBUST-CUTOFF / ROBUST-THRESH)
% and saves to `<savepath>/NN_lpcgc_comparison.fig`.
%
% T, D : cell arrays from run_lpcgcComparison (size ncutoff x nthresh).

    if ~exist('figure_number','var') || isempty(figure_number) || ~isnumeric(figure_number)
        figure_number = 0;
    end
    figure_number = figure_number + 1;

    params = analysis.params;
    ncutoff = length(params.lowpass_cutoffs);
    nthresh = length(params.saccade_threshs);

    ntColor = getOrDefault(params, 'ntColor', [1 0 0]);
    tnColor = getOrDefault(params, 'tnColor', [0 0 1]);
    exp_name = getOrDefault(params, 'exp_name', '');

    fig = uifigure('Units','normalized', 'Position',[0 0 1 1], 'WindowState','maximized');
    tg  = uitabgroup(fig, 'Units','normalized', 'Position',[0 0 1 1]);

    tb(1) = uitab(tg, 'Title','MAIN');
    h(1)  = tiledlayout(tb(1), ncutoff, nthresh, 'TileSpacing','tight', 'Padding','compact');
    title(h(1), sprintf("(%s) LOWPASS FREQUENCY VS SACCADE THRESHOLD", exp_name), 'FontSize',16, 'Interpreter','none');

    tb(2) = uitab(tg, 'Title','FILTER');
    h(2)  = tiledlayout(tb(2), ncutoff, nthresh, 'TileSpacing','tight', 'Padding','compact');
    title(h(2), sprintf("(%s) LOWPASS FREQUENCY COMPARISON BY SACCADE THRESHOLD", exp_name), 'FontSize',16, 'Interpreter','none');

    tb(3) = uitab(tg, 'Title','THRESH');
    h(3)  = tiledlayout(tb(3), ncutoff, nthresh, 'TileSpacing','tight', 'Padding','compact');
    title(h(3), sprintf("(%s) SACCADE THRESH COMPARISON BY LOWPASS FREQUENCY", exp_name), 'FontSize',16, 'Interpreter','none');

    snum = 0;
    S = struct;
    t = 1:length(D{1,1}(1).pooled_good_cyclemean_diff);
    ystim = mean([T{1,1}(1).pooled_drumvel_cyclemean_cvd.cycleMean; ...
                  T{1,1}(1).pooled_drumvel_cyclemean_cvd.cycleMean], 1, 'omitnan');

    % ---------------- TAB 1: MAIN ----------------
    axnum = 0;
    ax = gobjects(0);
    for ii = 1:ncutoff
        fcii = params.lowpass_cutoffs(ii);
        for jj = 1:nthresh
            thjj = params.saccade_threshs(jj);
            tpi = T{ii,jj};
            ddi_all = D{ii,jj};
            % Pair = (first timepoint) -> (last timepoint). Pick the matching
            % diffdata entry so post ~= last and the diff aligns.
            preIdx  = 1;
            postIdx = length(tpi);
            ppMat = vertcat(ddi_all.timePointPrePostIds);
            ddi   = ddi_all(ppMat(:,1) == preIdx & ppMat(:,2) == postIdx);

            ypre_cvd  = tpi(preIdx).pooled_good_cyclemean_cvd;
            ypre_cm   = tpi(preIdx).pooled_good_cyclemean_cm;
            ypre      = ypre_cvd.cycleMean;

            ypost_cvd = tpi(postIdx).pooled_good_cyclemean_cvd;
            ypost_cm  = tpi(postIdx).pooled_good_cyclemean_cm;
            ypost     = ypost_cvd.cycleMean;

            ydiff_cm  = ddi.pooled_good_cyclemean_diff_cm;
            ydiff     = ddi.pooled_good_cyclemean_diff;
            ydiff_sem = ddi.pooled_good_cyclemean_diff_SEM;

            % Pointwise two-sample Welch t-interval for ydiff
            nu = (ypre_cvd.cycleSEM.^2 + ypost_cvd.cycleSEM.^2).^2 ./ ...
                (ypre_cvd.cycleSEM.^4 ./ (ypre_cvd.nValidCyclesAt - 1) + ...
                 ypost_cvd.cycleSEM.^4 ./ (ypost_cvd.nValidCyclesAt - 1));
            alpha = ypre_cvd.alpha;
            tCrit = tinv(1 - alpha/2, nu);
            diff_ci = [ydiff - tCrit .* ydiff_sem; ydiff + tCrit .* ydiff_sem];

            axnum = axnum + 1;
            ax(axnum) = nexttile(h(1));
            maybeCustomToolbar(ax(axnum));
            hold(ax(axnum), 'on');
            yline(ax(axnum), 0, 'Color',[0 0 0]+0.7, 'LineWidth',0.1, 'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
            fill(ax(axnum), [t fliplr(t)], [diff_ci(1,:) fliplr(diff_ci(2,:))], ...
                'g', 'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
            plot(ax(axnum), t, ystim, '--k', 'LineWidth',2, 'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
            plot(ax(axnum), t, ypre,  '-b', 'LineWidth',2, 'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
            plot(ax(axnum), t, ypost, '-r', 'LineWidth',2, 'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
            plot(ax(axnum), t, ydiff, '-g', 'LineWidth',2, 'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
            hold(ax(axnum), 'off');
            xlim(ax(axnum), [t(1), t(min(501,end))]);
            if jj==1,       ylabel(ax(axnum), sprintf('%g HZ LOWPASS FILTERED', fcii)); end
            if ii==ncutoff, xlabel(ax(axnum), sprintf('%g SACCADE THRESHOLD', thjj)); end
            box(ax(axnum), 'on');

            S = appendTraceEntry(S, 'ypre',  fcii, thjj, ypre,  ypre_cvd.nCycles, [0 0 1], ypre_cm);
            snum = snum + 1;
            S = appendTraceEntry(S, 'ypost', fcii, thjj, ypost, ypost_cvd.nCycles, [1 0 0], ypost_cm);
            snum = snum + 1;
            S = appendTraceEntry(S, 'ydiff', fcii, thjj, ydiff, NaN,                 [0 1 0], ydiff_cm);
            snum = snum + 1;
        end
    end
    linkaxes(ax, 'xy');

    % ---------------- TAB 2: FILTER ----------------
    coltypes  = {'ypre', 'ypost', 'ydiff'};
    linetypes = {'-','--',':'};
    mkrtypes  = {'square','diamond','o'};
    axnum = 0;
    ax2 = gobjects(0);
    for jj = 1:nthresh
        thjj = num2str(params.saccade_threshs(jj));
        for ii = 1:length(coltypes)
            rtii = coltypes{ii};
            Sii  = S(strcmpi({S.type},rtii) & strcmpi({S.thresh},thjj));
            clrs = adjacentColors(Sii(1).color, 0.4);
            freqs = {Sii.freq};

            axnum = axnum + 1;
            ax2(axnum) = nexttile(h(2));
            maybeCustomToolbar(ax2(axnum));
            hold(ax2(axnum), 'on');
            yline(ax2(axnum), 0, 'Color',[0 0 0]+0.7, 'LineWidth',0.1, 'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
            plot(ax2(axnum), t, ystim, '--k', 'LineWidth',2, 'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
            for kk = 1:length(freqs)
                meanpt1 = Sii(kk).meanPeakTime1;
                medpt1  = Sii(kk).medPeakTime1;
                medpt2  = Sii(kk).medPeakTime2 + 500;
                xline(ax2(axnum), [medpt1, medpt2], 'k', '', 'LineStyle',':', 'LineWidth',1, 'Color',clrs{kk}, 'HandleVisibility','off');
                xline(ax2(axnum), [meanpt1, medpt2], 'k', '', 'LineStyle','--', 'LineWidth',1, 'Color',clrs{kk}, 'HandleVisibility','off');
                dntxt = sprintf('%s (%s Hz) | N = %s', rtii, Sii(kk).freq, Sii(kk).ncycles);
                pt1 = Sii(kk).peakTime1;
                plot(ax2(axnum), t, Sii(kk).data, linetypes{kk}, 'Color',clrs{kk}, 'LineWidth',2, 'DisplayName',dntxt, 'HitTest','off', 'PickableParts','none');
                plot(ax2(axnum), t(pt1), Sii(kk).data(pt1), mkrtypes{kk}, 'Color',clrs{kk}, 'MarkerFaceColor',clrs{kk}, 'MarkerSize',12, 'LineWidth',1, 'DisplayName',dntxt, 'HitTest','off', 'PickableParts','none');
            end
            hold(ax2(axnum), 'off');
            xlim(ax2(axnum), [t(1), t(min(501,end))]);
            if ii==1, ylabel(ax2(axnum), sprintf('%s SACCADE THRESHOLD', thjj)); end
            if jj==length(coltypes), xlabel(ax2(axnum), sprintf('%s', rtii)); end
            legend(ax2(axnum));
            box(ax2(axnum), 'on');
        end
    end
    linkaxes(ax2, 'xy');

    % ---------------- TAB 3: THRESH ----------------
    axnum = 0;
    ax3 = gobjects(0);
    for jj = 1:ncutoff
        fcjj = num2str(params.lowpass_cutoffs(jj));
        for ii = 1:length(coltypes)
            rtii = coltypes{ii};
            Sii = S(strcmpi({S.type},rtii) & strcmpi({S.freq},fcjj));
            clrs = adjacentColors(Sii(1).color, 0.4);
            threshs = {Sii.thresh};

            axnum = axnum + 1;
            ax3(axnum) = nexttile(h(3));
            maybeCustomToolbar(ax3(axnum));
            hold(ax3(axnum), 'on');
            yline(ax3(axnum), 0, 'Color',[0 0 0]+0.7, 'LineWidth',0.1, 'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
            plot(ax3(axnum), t, ystim, '--k', 'LineWidth',2, 'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
            for kk = 1:length(threshs)
                meanpt1 = Sii(kk).meanPeakTime1;
                medpt1  = Sii(kk).medPeakTime1;
                medpt2  = Sii(kk).medPeakTime2 + 500;
                xline(ax3(axnum), [medpt1, medpt2], 'k', '', 'LineStyle',':', 'LineWidth',1, 'Color',clrs{kk}, 'HandleVisibility','off');
                xline(ax3(axnum), [meanpt1, medpt2], 'k', '', 'LineStyle','--', 'LineWidth',1, 'Color',clrs{kk}, 'HandleVisibility','off');
                dntxt = sprintf('%s (thresh = %s) | N = %s', rtii, threshs{kk}, Sii(kk).ncycles);
                pt1 = Sii(kk).peakTime1;
                plot(ax3(axnum), t, Sii(kk).data, linetypes{kk}, 'Color',clrs{kk}, 'LineWidth',2, 'DisplayName',dntxt, 'HitTest','off', 'PickableParts','none');
                plot(ax3(axnum), t(pt1), Sii(kk).data(pt1), mkrtypes{kk}, 'Color',clrs{kk}, 'MarkerFaceColor',clrs{kk}, 'MarkerSize',12, 'LineWidth',1, 'DisplayName',dntxt, 'HitTest','off', 'PickableParts','none');
            end
            hold(ax3(axnum), 'off');
            xlim(ax3(axnum), [t(1), t(min(501,end))]);
            if ii==1, ylabel(ax3(axnum), sprintf('%s LOWPASS FREQUENCY', fcjj)); end
            if jj==length(coltypes), xlabel(ax3(axnum), sprintf('%s', rtii)); end
            legend(ax3(axnum));
            box(ax3(axnum), 'on');
        end
    end
    linkaxes(ax3, 'xy');

    % ---------------- TAB 4: ROBUST-CUTOFF ----------------
    tb(4) = uitab(tg, 'Title','ROBUST-CUTOFF');
    h(4)  = tiledlayout(tb(4), 2, 2, 'TileSpacing','tight', 'Padding','compact');
    title(h(4), sprintf("(%s) ROBUSTNESS-TO-CUTOFF BY SACCADE THRESHOLD", exp_name), 'FontSize',16, 'Interpreter','none');

    mRC = compute_robustness_metrics(S, params.saccade_threshs, 'thresh', 'freq');
    ax_rc(1) = nexttile(h(4));
    plot_robustness_subplot(ax_rc(1), params.saccade_threshs, mRC.argmax_NT, mRC.argmax_TN, 'argmax spread (ms)', 'saccade threshold', true, ntColor, tnColor);
    legend(ax_rc(1), 'Location','best');
    ax_rc(2) = nexttile(h(4));
    plot_robustness_subplot(ax_rc(2), params.saccade_threshs, mRC.cmean_NT, mRC.cmean_TN, 'cmean spread (ms)', 'saccade threshold', true, ntColor, tnColor);
    ax_rc(3) = nexttile(h(4));
    plot_robustness_subplot(ax_rc(3), params.saccade_threshs, mRC.cmed_NT, mRC.cmed_TN, 'cmed spread (ms)', 'saccade threshold', true, ntColor, tnColor);
    ax_rc(4) = nexttile(h(4));
    plot_robustness_subplot(ax_rc(4), params.saccade_threshs, mRC.auc, [], 'band-AUC (deg/s \cdot ms)', 'saccade threshold', false, ntColor, tnColor);

    % ---------------- TAB 5: ROBUST-THRESH ----------------
    tb(5) = uitab(tg, 'Title','ROBUST-THRESH');
    h(5)  = tiledlayout(tb(5), 2, 2, 'TileSpacing','tight', 'Padding','compact');
    title(h(5), sprintf("(%s) ROBUSTNESS-TO-THRESHOLD BY LOWPASS FREQUENCY", exp_name), 'FontSize',16, 'Interpreter','none');

    mRT = compute_robustness_metrics(S, params.lowpass_cutoffs, 'freq', 'thresh');
    ax_rt(1) = nexttile(h(5));
    plot_robustness_subplot(ax_rt(1), params.lowpass_cutoffs, mRT.argmax_NT, mRT.argmax_TN, 'argmax spread (ms)', 'lowpass cutoff (Hz)', true, ntColor, tnColor);
    legend(ax_rt(1), 'Location','best');
    ax_rt(2) = nexttile(h(5));
    plot_robustness_subplot(ax_rt(2), params.lowpass_cutoffs, mRT.cmean_NT, mRT.cmean_TN, 'cmean spread (ms)', 'lowpass cutoff (Hz)', true, ntColor, tnColor);
    ax_rt(3) = nexttile(h(5));
    plot_robustness_subplot(ax_rt(3), params.lowpass_cutoffs, mRT.cmed_NT, mRT.cmed_TN, 'cmed spread (ms)', 'lowpass cutoff (Hz)', true, ntColor, tnColor);
    ax_rt(4) = nexttile(h(5));
    plot_robustness_subplot(ax_rt(4), params.lowpass_cutoffs, mRT.auc, [], 'band-AUC (deg/s \cdot ms)', 'lowpass cutoff (Hz)', false, ntColor, tnColor);

    % ---------------- Save ----------------
    if ~isequal(savepath,0) && isfolder(savepath)
        filename = sprintf('%02d_lpcgc_comparison.fig', figure_number);
        savefig(fig, fullfile(savepath, filename));
        fprintf('    Saved: %s\n', filename);
    end
end


%% ---------------- Local helpers ----------------

function S = appendTraceEntry(S, typ, fcii, thjj, traceData, nCycles, color, cm)
    idx = numel(S) + 1;
    if numel(S) == 1 && ~isfield(S, 'type'), idx = 1; end
    S(idx).type = typ;
    S(idx).freq = num2str(fcii);
    S(idx).thresh = num2str(thjj);
    S(idx).data = traceData;
    if isnan(nCycles), S(idx).ncycles = 'NA'; else, S(idx).ncycles = num2str(nCycles); end
    S(idx).color = color;
    [~,pt1] = max(traceData(1:min(500,end)));
    [~,pt2] = min(traceData(min(501,end):end));
    S(idx).peakTime1      = pt1;
    S(idx).peakTime2      = pt2;
    S(idx).cmPeakTime1    = cm.eye.peak1TimeMs(1) + 1;
    S(idx).cmPeakTime2    = cm.eye.peak1TimeMs(2) + 1;
    S(idx).meanPeakTime1  = cm.eye.centroidMean(1);
    S(idx).meanPeakTime2  = cm.eye.centroidMean(2);
    S(idx).medPeakTime1   = cm.eye.centroidMedian(1);
    S(idx).medPeakTime2   = cm.eye.centroidMedian(2);
end

function cout = adjacentColors(base, offset)
    if nargin < 2, offset = 0.3; end
    idx = find(base == 0);
    c1 = base; c1(idx(1)) = offset;
    c2 = base; c2(idx(2)) = offset;
    cout = {c1, base, c2};
end

function m = compute_robustness_metrics(S, xVals, heldField, sweptField) %#ok<INUSD>
    traceTypes = {'ypre','ypost','ydiff'};
    nX = length(xVals);
    nT = length(traceTypes);
    m.argmax_NT = nan(nX, nT);  m.argmax_TN = nan(nX, nT);
    m.cmean_NT  = nan(nX, nT);  m.cmean_TN  = nan(nX, nT);
    m.cmed_NT   = nan(nX, nT);  m.cmed_TN   = nan(nX, nT);
    m.auc       = nan(nX, nT);
    for jj = 1:nX
        heldVal = num2str(xVals(jj));
        for tt = 1:nT
            rtii = traceTypes{tt};
            Sii = S(strcmpi({S.type},rtii) & strcmpi({S.(heldField)},heldVal));
            if numel(Sii) ~= 3
                warning('compute_robustness_metrics: expected 3 siblings, got %d (%s=%s, type=%s)', ...
                    numel(Sii), heldField, heldVal, rtii);
                continue;
            end
            m.argmax_NT(jj,tt) = max([Sii.peakTime1])     - min([Sii.peakTime1]);
            m.argmax_TN(jj,tt) = max([Sii.peakTime2])     - min([Sii.peakTime2]);
            m.cmean_NT(jj,tt)  = max([Sii.meanPeakTime1]) - min([Sii.meanPeakTime1]);
            m.cmean_TN(jj,tt)  = max([Sii.meanPeakTime2]) - min([Sii.meanPeakTime2]);
            m.cmed_NT(jj,tt)   = max([Sii.medPeakTime1])  - min([Sii.medPeakTime1]);
            m.cmed_TN(jj,tt)   = max([Sii.medPeakTime2])  - min([Sii.medPeakTime2]);
            traces = [Sii(1).data; Sii(2).data; Sii(3).data];
            m.auc(jj,tt) = sum(max(traces,[],1) - min(traces,[],1));
        end
    end
end

function plot_robustness_subplot(ax, xVals, dataNT, dataTN, yLab, xLab, isPeakTime, ntColor, tnColor)
    hold(ax, 'on');
    nX = length(xVals);
    catPos = 1:(nX+1);
    xOffset = 0.08;
    if isPeakTime
        plot_mean_ci_cats(ax, catPos - xOffset, dataNT, ntColor, 'NT');
        plot_mean_ci_cats(ax, catPos + xOffset, dataTN, tnColor, 'TN');
    else
        plot_mean_ci_cats(ax, catPos, dataNT, ntColor, 'band-AUC');
    end
    hold(ax, 'off');
    xlabel(ax, xLab); ylabel(ax, yLab);
    xticks(ax, catPos);
    tickLabels = [arrayfun(@(v) num2str(v), xVals(:)', 'UniformOutput', false), {'mean'}];
    xticklabels(ax, tickLabels);
    xlim(ax, [catPos(1)-0.5, catPos(end)+0.5]);
    xline(ax, catPos(end)-0.5, ':', 'Color',[0.5 0.5 0.5], 'HandleVisibility','off');
    box(ax, 'on'); grid(ax, 'on');
end

function plot_mean_ci_cats(ax, xPos, data, color, displayName)
    nX = size(data, 1);
    mu = nan(nX+1, 1);
    ci = nan(nX+1, 1);
    for jj = 1:nX
        vals = data(jj,:);
        [mu(jj), ci(jj)] = mean_with_ci(vals(~isnan(vals)));
    end
    flat = data(:);
    [mu(nX+1), ci(nX+1)] = mean_with_ci(flat(~isnan(flat)));
    plot(ax, xPos(1:nX), mu(1:nX), '-', 'Color',color, 'LineWidth',1.25, 'HandleVisibility','off');
    errorbar(ax, xPos, mu, ci, 'o', ...
        'MarkerFaceColor',color, 'MarkerEdgeColor',color, ...
        'Color',color, 'LineWidth',1, 'MarkerSize',7, 'CapSize',8, ...
        'LineStyle','none', 'DisplayName',displayName);
end

function [mu, ci] = mean_with_ci(vals)
    n = numel(vals);
    if n < 2, mu = mean(vals); ci = 0; return; end
    mu  = mean(vals);
    sem = std(vals) ./ sqrt(n);
    ci  = tinv(0.975, n-1) .* sem;
end

function v = getOrDefault(s, field, def)
    if isfield(s, field) && ~isempty(s.(field)), v = s.(field); else, v = def; end
end

function maybeCustomToolbar(ax)
    if exist('axesCustomToolbarButtons', 'file') == 2
        try, axesCustomToolbarButtons(ax); catch, end
    end
end
