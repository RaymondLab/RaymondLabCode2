clear; clc; close all;
delete(timerfindall);  % Clean up any orphaned timers from previous runs
params.script_folderpath = fileparts(which(matlab.desktop.editor.getActiveFilename));
addpath(genpath(params.script_folderpath));

% Note: This script will not work if the single_NDD_okrAnalysis was set up
% with "params.timepoint_groups" containing more than two groups (0 and 60)
save_filename = 'group_NDD_analysis.xlsx';

subconds = {'NR', '23DMR', '24DMR', 'DR'};
groups_to_plot = {'TNOR', 'TDIM', 'TDMR', '2TNOR', '2TDIM', '2TDMR'};
group_colors = {'-r', '-g', '-b', '-m', '-y', '-c'};
group_colors2 = {[0.7 0 0], [0 0.7 0], [0 0 0.7], [0.7 0 0.7], [0.7 0.7 0], [0 0.7 0.7]};
group_colors3 = {[0.5 0 0], [0 0.5 0], [0 0 0.5], [0.5 0 0.5], [0.5 0.5 0], [0 0.5 0.5]};
cond_colors = {'-r', '-g', '-m', '-b'};
cond_colors2 = {[0.7 0 0], [0 0.7 0], [0.7 0 0.7], [0 0 0.7]};
cond_colors3 = {[0.5 0 0], [0 0.5 0], [0.5 0 0.5], [0 0 0.5]};
jcolors = {'-b', '-r', '-g'};
jcolors2 = {[0 0 0.5], [0.5 0 0], [0 0.5 0]};
fs = 1000;
freq = 1;
cycleLength = round(fs / freq);
halfCycle = round(cycleLength / 2);
fignum = 0;


%% SCRIPT START
fprintf('\n--- SCRIPT START! ---\n');

% Get the folder to retrieve group data from
disp('SELECT ROOT FOLDER THAT CONTAINS ALL ANALYZED DATA...');
source_folder = uigetdir(pwd, 'Select folder containing group data');
if isequal(source_folder, 0)
    error('No source folder selected. Aborting script!');
elseif ~isfolder(source_folder)
    error('Source folder does not exist: %s', source_folder);
end
fprintf('\n   Source data folder path selected: %s\n', source_folder);

% Get the folder to save group data to
disp('SELECT FOLDER TO SAVE ANALYSIS RESULTS TO...');
save_folderpath = uigetdir(source_folder, 'Select folder to save analysis to');
if isequal(save_folderpath, 0)
    error('No save folder selected. Aborting script!');
elseif ~isfolder(save_folderpath)
    error('Save folder does not exist: %s', save_folderpath);
end
fprintf('\n   Save folder path selected: %s\n', save_folderpath);

% Get experiment information for all analysis results in source folder
[expinfo,nexp] = get_expInfo(source_folder);

disp('Retrieving analysis data from experiments...');
% Retrieve data from each experiment file
xlsx = struct;  % Preallocate struct for storing tables for excel sheet
for i = 1:nexp
    fp = fullfile(expinfo(i).folder, expinfo(i).name);
    DD = readtable(fp, 'Sheet','DiffData');
    expinfo(i).metrics = readtable(fp, 'Sheet','DiffInfo');
    fields = DD.Properties.VariableNames;
    for j = 1:numel(fields)
        expinfo(i).(fields{j}) = DD.(fields{j});
        if ~ismember(fields{j}, {'CycleTimes','StimVel'})
            if ~isfield(xlsx, fields{j})
                xlsx.(fields{j}) = table;
                if ~contains(fields{j}, 'SEM')
                    xlsx.(fields{j}).StimVel = DD.StimVel(:);
                end
            end
        end
    end
end

disp('Compiling data for excel sheet...');
% Compile struct for the output table of results
fields = fieldnames(xlsx);
for j = 1:nexp
    for i = 1:numel(fields)
        ifield = fields{i};
        jheader = strrep(expinfo(j).name, '_analysis-diffdata.xlsx', '');
        xlsx.(ifield).(jheader) = expinfo(j).(ifield);
    end
end

% Save results to excel file at the chosen destination
if ~isfolder(save_folderpath)
    error('Chosen destination folder is invalid: %s', save_folderpath);
end

disp('Saving data to excel sheet at chosen destination...');
xlsx_filepath = fullfile(save_folderpath, save_filename);
if isfile(xlsx_filepath), delete(xlsx_filepath); end
for i = 1:numel(fields)
    ifield = fields{i};
    writetable(xlsx.(ifield), xlsx_filepath, 'Sheet',ifield);
end

disp('Plotting group analysis results...');

% Get unique condition types across all experiments
uniqueConds = intersect(subconds, {expinfo.cond}, 'stable');
nConds = numel(uniqueConds);
fprintf('Found %d unique condition(s): %s\n', nConds, strjoin(uniqueConds, ', '));

% Create single summary figure with tabs for each condition
fignum = fignum + 1;
summaryFigIdx = fignum;
fig(fignum) = uifigure('Units','normalized', 'Position',[0 0 1 1], 'WindowState','maximized');
summaryTg = uitabgroup(fig(fignum), 'Units','normalized', 'Position',[0 0 1 1]);
summary_fig_filename = sprintf('%02d_Summary_TimePoint_groupAnalysis.fig', fignum);

% Preallocate global structs for storing group arrays across all conditions
all_group_cyclemeans = struct;
allGcmIdx = 0;
allGcmUsedKeys = {};
allGdIdx = 0;
all_grpDiffs = struct;

% =========================================================================
% Loop over each unique mouse condition
% =========================================================================
ngroups = length(groups_to_plot);
sp_per_row = 3;
for ic = 1:nConds
    icond = uniqueConds{ic};
    fprintf('\n--- Processing mouse condition: %s ---\n', icond);

    % Filter experiments to this condition
    condExp = expinfo(strcmpi({expinfo.cond}, icond));

    % Get all the timepoint types from experiments in this condition
    alltypes = get_orderedRowTypes(condExp);

    % Plot comparison of timepoints in separate figures per group
    assert(mod(numel(alltypes), 3) == 0, ...
        'alltypes length (%d) is not divisible by 3 for condition "%s".', numel(alltypes), icond);
    ntpgroups = round(numel(alltypes) / 3);
    alltypes2 = reshape(alltypes, [3, ntpgroups]);
    axs = gobjects(0);  % Reset axes array for this condition
    for i = 1:ngroups
        igroup = groups_to_plot{i};
        iexp = condExp(strcmpi({condExp.task}, igroup));
        nGroupExp = length(iexp);

        % Skip groups with no matching experiments
        if nGroupExp == 0
            warning('No experiments found for group "%s" in condition "%s". Skipping.', igroup, icond);
            continue;
        end

        % Set up plot
        fignum = fignum + 1;
        fig(fignum) = uifigure('Units','normalized', 'Position',[0 0 1 1], 'WindowState','maximized');
        tg = uitabgroup(fig(fignum), 'Units','normalized', 'Position',[0 0 1 1]);
        fig_filename = sprintf('%02d_%s_%s_exps_groupAnalysis.fig', fignum, icond, igroup);

        for tpg = 1:ntpgroups
            tps = alltypes2(:,tpg);
            ntps = length(tps);
            
            tb(tpg) = uitab(tg, "Title",tps{end});
            h = tiledlayout(tb(tpg), 2, 2, 'TileSpacing','compact', 'Padding','compact');
            title(h, sprintf('[%s] %s Experiments Analysis', icond, igroup), 'FontSize',16, 'Interpreter','none');
            
            % Compute once per tab group, not per timepoint
            cycleTimes = iexp(1).CycleTimes;
            group_StimVel = mean([iexp.StimVel], 2, 'omitnan');

            tp_cyclemeans = zeros(cycleLength, ntps);
            tp_res = cell(ntps, 1);
            for j = 1:ntps
                jtps = tps{j};
                jcyclemeans = zeros(cycleLength, nGroupExp);
        
                axs(i,j) = nexttile(h);
                axesCustomToolbarButtons(axs(i,j));
                hold(axs(i,j), 'on');
                yline(axs(i,j), 0, 'Color',[0 0 0]+0.7, 'LineWidth',0.1, 'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
                plot(axs(i,j), cycleTimes, group_StimVel, '--k', 'LineWidth',2, 'DisplayName','Stimulus Velocity', 'HitTest','off', 'PickableParts','none');
                
                for k = 1:nGroupExp
                    jcyclemeans(:,k) = iexp(k).(jtps);
                    plot(axs(i,j), cycleTimes, jcyclemeans(:,k), 'LineWidth',1, 'DisplayName',iexp(k).sub, 'HitTest','off', 'PickableParts','none');
                end
                tp_cyclemeans(:,j) = mean(jcyclemeans,2,'omitnan');  % Group cycle-mean
                
                resjj = calc_cycleMetrics(tp_cyclemeans(:,j), group_StimVel);
                rj = parse_cycleMetricsResults(resjj, halfCycle, fs);
        
                xline(axs(i,j), rj.meanCenAdj, 'k', rj.meanTxt, 'LineStyle',':', 'LineWidth',0.001, 'LabelVerticalAlignment','middle', 'LabelOrientation','horizontal', 'Color',jcolors2{j}, 'HandleVisibility','off');
                xline(axs(i,j), rj.medCenAdj, 'k', rj.medTxt, 'LineStyle','--', 'LineWidth',0.001, 'LabelVerticalAlignment','bottom', 'LabelOrientation','horizontal', 'Color',jcolors2{j}, 'HandleVisibility','off');
                xline(axs(i,j), rj.meanCenAdj, jcolors{j}, 'LineStyle',':', 'LineWidth',1.5, 'HandleVisibility','off');
                xline(axs(i,j), rj.medCenAdj, jcolors{j}, 'LineStyle','--', 'LineWidth',1, 'HandleVisibility','off');
                plot(axs(i,j), cycleTimes, tp_cyclemeans(:,j), '-k', 'LineWidth',3, 'DisplayName','Mean', 'HitTest','off', 'PickableParts','none');
                plot(axs(i,j), rj.ntLocs, rj.ntPks, 'm.', 'MarkerSize',20, 'DisplayName','Peaks');
                plot(axs(i,j), rj.tnLocs, rj.tnPks, 'm.', 'MarkerSize',20, 'HandleVisibility','off');
                text(axs(i,j), rj.ntLocs, rj.ntPks, cellstr(" "+string(rj.ntLocsMs))+"ms", 'Color','m', 'VerticalAlignment','bottom');
                text(axs(i,j), rj.tnLocs, rj.tnPks, cellstr(" "+string(rj.tnLocsMs))+"ms", 'Color','m', 'VerticalAlignment','top');
                hold(axs(i,j), 'off');
                xlim(axs(i,j), [0 cycleTimes(end)]);
                title(axs(i,j), sprintf('%s', jtps), 'FontSize',12, 'Interpreter','none');
                legend(axs(i,j));
                box(axs(i,j), 'on');
        
                tp_res{j} = rj;
            end
        
            lastSP = j+1;
            axs(i,lastSP) = nexttile(h);
            axesCustomToolbarButtons(axs(i,lastSP));
            hold(axs(i,lastSP), 'on');
            yline(axs(i,lastSP), 0, 'Color',[0 0 0]+0.7, 'LineWidth',0.1, 'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
            plot(axs(i,lastSP), cycleTimes, group_StimVel, '--k', 'LineWidth',2, 'DisplayName','Stimulus Velocity', 'HitTest','off', 'PickableParts','none');
            for j = 1:ntps
                jtps = tps{j};
                meanIds = max(1, min(cycleLength, round(tp_res{j}.meanCenAdj * fs) + 1));
                medIds = max(1, min(cycleLength, round(tp_res{j}.medCenAdj * fs) + 1));
                jdname = sprintf('%s Timepoint Group Cycle-mean', jtps);
                plot(axs(i,lastSP), cycleTimes, tp_cyclemeans(:,j), jcolors{j}, 'LineWidth',2, 'DisplayName',jdname, 'HitTest','off', 'PickableParts','none');
                plot(axs(i,lastSP), cycleTimes(meanIds), tp_cyclemeans(meanIds,j), '^', 'Color',jcolors2{j}, 'MarkerFaceColor',jcolors2{j}, 'MarkerSize',8, 'HandleVisibility','off');
                plot(axs(i,lastSP), cycleTimes(medIds), tp_cyclemeans(medIds,j), 'o', 'Color',jcolors2{j}, 'MarkerSize',10, 'HandleVisibility','off');
            end
            amv = max(abs(group_StimVel));

            % Guard against fewer than 2 timepoints
            if ntps >= 2
                meanDiffsSecs = tp_res{2}.meanCenAdj - tp_res{1}.meanCenAdj;
                meanDiffs = round(meanDiffsSecs * fs);
                meanHalfDiffs = abs(meanDiffsSecs) ./ 2;
                mean_xmid = (tp_res{1}.meanCenAdj + tp_res{2}.meanCenAdj) ./ 2;
                mean_ymid = [amv+(0.2*amv), (-amv)-(0.2*amv)];
                errorbar(axs(i,lastSP), mean_xmid, mean_ymid, meanHalfDiffs, 'horizontal', 'Color','k', 'LineWidth',1, 'LineStyle','none', 'HandleVisibility','off');
                text(axs(i,lastSP), mean_xmid+meanHalfDiffs, mean_ymid, ...
                    cellstr("    Mean\Deltat: "+string(meanDiffs))+"ms", 'FontSize',12, 'HorizontalAlignment','left', 'VerticalAlignment','middle');
            
                medDiffsSecs = tp_res{2}.medCenAdj - tp_res{1}.medCenAdj;
                medDiffs = round(medDiffsSecs * fs);
                medHalfDiffs = abs(medDiffsSecs) ./ 2;
                med_xmid = (tp_res{1}.medCenAdj + tp_res{2}.medCenAdj) ./ 2;
                med_ymid = [amv+(0.4*amv), (-amv)-(0.4*amv)];
                errorbar(axs(i,lastSP), med_xmid, med_ymid, medHalfDiffs, 'horizontal', 'Color','k', 'LineWidth',1, 'LineStyle','none', 'HandleVisibility','off');
                text(axs(i,lastSP), med_xmid+medHalfDiffs, med_ymid, ...
                    cellstr("    Median\Deltat: "+string(medDiffs))+"ms", 'FontSize',12, 'HorizontalAlignment','left', 'VerticalAlignment','middle');
            end
        
            hold(axs(i,lastSP), 'off');
            xlim(axs(i,lastSP), [0 cycleTimes(end)]);
            title(axs(i,lastSP), sprintf('%s Timepoints Cycle-Mean Comparison', igroup), 'FontSize',12, 'Interpreter','none');
            legend(axs(i,lastSP), 'Interpreter','none');
            box(axs(i,lastSP), 'on');
        
            % Add results to the global group struct
            for j = 1:ntps
                jtps = tps{j};
                pairKey = strcat(icond, '||', igroup, '||', jtps);
                if ismember(pairKey, allGcmUsedKeys)
                    continue;
                end
                allGcmUsedKeys{end+1} = pairKey; %#ok<AGROW>
                allGcmIdx = allGcmIdx + 1;
                all_group_cyclemeans(allGcmIdx).cond = icond;
                all_group_cyclemeans(allGcmIdx).group = igroup;
                all_group_cyclemeans(allGcmIdx).types = jtps;
                all_group_cyclemeans(allGcmIdx).cycletimes = cycleTimes;
                all_group_cyclemeans(allGcmIdx).stimvels = group_StimVel;
                all_group_cyclemeans(allGcmIdx).cyclemeans = tp_cyclemeans(:,j);
                all_group_cyclemeans(allGcmIdx).cyclemetrics = tp_res{j};
            end

            if ntps >= 2
                allGdIdx = allGdIdx + 1;
                all_grpDiffs(allGdIdx).cond = icond;
                all_grpDiffs(allGdIdx).group = igroup;
                all_grpDiffs(allGdIdx).dt_type = tps{end};
                all_grpDiffs(allGdIdx).dt_mean = meanDiffs(:);
                all_grpDiffs(allGdIdx).dt_median = medDiffs(:);
            end
        end
        linkaxes(axs(i,:));
        drawnow;

        if ~isequal(fig_filename, 0) && isfolder(save_folderpath)
            savefig(fig(fignum), fullfile(save_folderpath, fig_filename));
            fprintf('    Figure %d saved as: %s\n', fignum, fig_filename);
        end
    end

    % After group loop - get global min/max across all subplots for this condition
    allYLims = vertcat(axs.YLim);  % n x 2 matrix
    globalYMin = min(allYLims(:,1));
    globalYMax = max(allYLims(:,2));
    globalYPad = min([abs(globalYMin), abs(globalYMax)]);

    % Filter global structs to this condition for the summary tab
    condMask = strcmpi({all_group_cyclemeans.cond}, icond);
    cond_gcm = all_group_cyclemeans(condMask);
    condDiffMask = strcmpi({all_grpDiffs.cond}, icond);
    cond_grpDiffs = all_grpDiffs(condDiffMask);

    % Plot comparisons of groups by timepoint for this condition
    cycleTimes = cond_gcm(1).cycletimes;
    grandStimVel = mean([cond_gcm.stimvels], 2, 'omitnan');

    % Update alltypes array to only unique types to avoid redundancies
    alltypes = get_orderedRowTypesV2(condExp);
    ntypes = numel(alltypes);

    summaryTb = uitab(summaryTg, 'Title', icond);
    h2 = tiledlayout(summaryTb, ceil(ntypes/sp_per_row), sp_per_row, 'TileSpacing','compact', 'Padding','compact');
    title(h2, sprintf('[%s] Group Comparisons by Timepoints', icond), 'FontSize',16, 'Interpreter','none');

    ax = gobjects(1, ntypes);
    for i = 1:ntypes
        itype = alltypes{i};
        igcms = cond_gcm(ismember({cond_gcm.types}, itype));

        ax(i) = nexttile(h2);
        axesCustomToolbarButtons(ax(i));
        hold(ax(i), 'on');
        yline(ax(i), 0, 'Color',[0 0 0]+0.7, 'LineWidth',0.1, 'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
        plot(ax(i), cycleTimes, grandStimVel, '--k', 'LineWidth',2, 'DisplayName','Stimulus Velocity', 'HitTest','off', 'PickableParts','none');

        ngrps = length(igcms);
        for j = 1:ngrps
            jgrp = igcms(j).group;
            jColorIdx = find(strcmpi(groups_to_plot, jgrp), 1);
            meanIds = max(1, min(cycleLength, round(igcms(j).cyclemetrics.meanCenAdj * fs) + 1));
            medIds = max(1, min(cycleLength, round(igcms(j).cyclemetrics.medCenAdj * fs) + 1));
            plot(ax(i), cycleTimes, igcms(j).cyclemeans, group_colors{jColorIdx}, 'LineWidth',2, 'DisplayName',jgrp, 'HitTest','off', 'PickableParts','none');
            plot(ax(i), cycleTimes(meanIds), igcms(j).cyclemeans(meanIds), '^', 'Color',group_colors2{jColorIdx}, 'MarkerFaceColor',group_colors2{jColorIdx}, 'MarkerSize',8, 'HandleVisibility','off');
            plot(ax(i), cycleTimes(medIds), igcms(j).cyclemeans(medIds), 'o', 'Color',group_colors2{jColorIdx}, 'MarkerSize',10, 'HandleVisibility','off');
            if contains(itype, 'Diff_')
                ijgd = cond_grpDiffs(ismember({cond_grpDiffs.group},jgrp) & contains({cond_grpDiffs.dt_type},itype));
                if ~isempty(ijgd)
                    ijMeanDiff = ijgd.dt_mean;
                    ijMedDiff = ijgd.dt_median;
                    ySpacing = 0.2 * j;
                    ytxt = [0 - (ySpacing*globalYPad), 0 + ((0.75 - ySpacing)*globalYPad)];
                    nl = newline;
                    txtstrs = cellstr(string(jgrp) + nl + "    Mean\Deltat: " + string(ijMeanDiff) + "ms" + " | Median\Deltat: " + string(ijMedDiff) + "ms");
                    text(ax(i), cycleTimes([round(0.05*cycleLength),round(0.55*cycleLength)]), ytxt, ...
                        txtstrs, 'Color',group_colors3{jColorIdx}, 'FontSize',10, 'HorizontalAlignment','left', 'VerticalAlignment','middle');
                end
            end
        end

        hold(ax(i), 'off');
        xlim(ax(i), [0 cycleTimes(end)]);
        title(ax(i), sprintf('%s Timepoints by Group', itype), 'FontSize',12, 'Interpreter','none');
        legend(ax(i), 'Interpreter','none');
        box(ax(i), 'on');
    end
    linkaxes(ax);

    fprintf('--- Finished mouse condition: %s ---\n', icond);
end
% =========================================================================
% End of condition loop
% =========================================================================

% =========================================================================
% Loop over each task group — tabs showing conditions as series
% =========================================================================
fprintf('\n--- Generating per-task-group summary tabs ---\n');
for ig = 1:ngroups
    igroup = groups_to_plot{ig};

    % Filter global structs to this task group
    grpMask = strcmpi({all_group_cyclemeans.group}, igroup);
    if ~any(grpMask)
        warning('No data found for task group "%s". Skipping.', igroup);
        continue;
    end
    grp_gcm = all_group_cyclemeans(grpMask);
    grpDiffMask = strcmpi({all_grpDiffs.group}, igroup);
    grp_grpDiffs = all_grpDiffs(grpDiffMask);

    % Get ordered unique types for this group's data
    grpExp = expinfo(strcmpi({expinfo.task}, igroup));
    alltypes = get_orderedRowTypesV2(grpExp);
    ntypes = numel(alltypes);

    cycleTimes = grp_gcm(1).cycletimes;
    grandStimVel = mean([grp_gcm.stimvels], 2, 'omitnan');

    % Compute Y-axis padding from this group's data
    grpYLims = [min(cellfun(@min, {grp_gcm.cyclemeans})), max(cellfun(@max, {grp_gcm.cyclemeans}))];
    grpYPad = min(abs(grpYLims));

    summaryTb = uitab(summaryTg, 'Title', igroup);
    h2 = tiledlayout(summaryTb, ceil(ntypes/sp_per_row), sp_per_row, 'TileSpacing','compact', 'Padding','compact');
    title(h2, sprintf('[%s] Condition Comparisons by Timepoints', igroup), 'FontSize',16, 'Interpreter','none');

    ax2 = gobjects(1, ntypes);
    for i = 1:ntypes
        itype = alltypes{i};
        igcms = grp_gcm(ismember({grp_gcm.types}, itype));

        ax2(i) = nexttile(h2);
        axesCustomToolbarButtons(ax2(i));
        hold(ax2(i), 'on');
        yline(ax2(i), 0, 'Color',[0 0 0]+0.7, 'LineWidth',0.1, 'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
        plot(ax2(i), cycleTimes, grandStimVel, '--k', 'LineWidth',2, 'DisplayName','Stimulus Velocity', 'HitTest','off', 'PickableParts','none');

        nconds_here = length(igcms);
        for j = 1:nconds_here
            jcond = igcms(j).cond;
            jColorIdx = find(strcmpi(uniqueConds, jcond), 1);
            meanIds = max(1, min(cycleLength, round(igcms(j).cyclemetrics.meanCenAdj * fs) + 1));
            medIds = max(1, min(cycleLength, round(igcms(j).cyclemetrics.medCenAdj * fs) + 1));
            plot(ax2(i), cycleTimes, igcms(j).cyclemeans, cond_colors{jColorIdx}, 'LineWidth',2, 'DisplayName',jcond, 'HitTest','off', 'PickableParts','none');
            plot(ax2(i), cycleTimes(meanIds), igcms(j).cyclemeans(meanIds), '^', 'Color',cond_colors2{jColorIdx}, 'MarkerFaceColor',cond_colors2{jColorIdx}, 'MarkerSize',8, 'HandleVisibility','off');
            plot(ax2(i), cycleTimes(medIds), igcms(j).cyclemeans(medIds), 'o', 'Color',cond_colors2{jColorIdx}, 'MarkerSize',10, 'HandleVisibility','off');
            if contains(itype, 'Diff_')
                ijgd = grp_grpDiffs(ismember({grp_grpDiffs.cond},jcond) & contains({grp_grpDiffs.dt_type},itype));
                if ~isempty(ijgd)
                    ijMeanDiff = ijgd.dt_mean;
                    ijMedDiff = ijgd.dt_median;
                    ySpacing = 0.2 * j;
                    ytxt = [0 - (ySpacing*grpYPad), 0 + ((0.75 - ySpacing)*grpYPad)];
                    nl = newline;
                    txtstrs = cellstr(string(jcond) + nl + "    Mean\Deltat: " + string(ijMeanDiff) + "ms" + " | Median\Deltat: " + string(ijMedDiff) + "ms");
                    text(ax2(i), cycleTimes([round(0.05*cycleLength),round(0.55*cycleLength)]), ytxt, ...
                        txtstrs, 'Color',cond_colors3{jColorIdx}, 'FontSize',10, 'HorizontalAlignment','left', 'VerticalAlignment','middle');
                end
            end
        end

        hold(ax2(i), 'off');
        xlim(ax2(i), [0 cycleTimes(end)]);
        title(ax2(i), sprintf('%s Timepoints by Condition', itype), 'FontSize',12, 'Interpreter','none');
        legend(ax2(i), 'Interpreter','none');
        box(ax2(i), 'on');
    end
    linkaxes(ax2);

    fprintf('    Added summary tab for task group: %s\n', igroup);
end
% =========================================================================
% End of task group loop
% =========================================================================

% Save the summary figure
if isfolder(save_folderpath)
    savefig(fig(summaryFigIdx), fullfile(save_folderpath, summary_fig_filename));
    fprintf('    Summary figure saved as: %s\n', summary_fig_filename);
end

% Delay re-arrangement to ensure all figures are rendered; clean up timer
t = timer('StartDelay', 3, 'ExecutionMode','singleShot', ...
    'TimerFcn', @(tmr,~) bringFigsToFront(fig, tmr));
start(t);

fprintf('\n--- Script completed!---\n\n');


%% Helper Functions
function bringFigsToFront(fig, tmr)
    for ii = length(fig):-1:1
        if isa(fig(ii), 'matlab.ui.Figure')
            fig(ii).Visible = 'off';
            fig(ii).Visible = 'on';
        else
            figure(fig(ii));
        end
    end
    uialert(fig(1), 'Analysis complete!', 'Done');
    stop(tmr);
    delete(tmr);
end

function ort = get_orderedRowTypesV2(expinfo)
    rt = arrayfun(@(x) x.metrics.Type, expinfo, 'UniformOutput', false);
    rt = unique(vertcat(rt{:}), 'stable');

    % Separate T entries and Diff entries
    isDiff = startsWith(rt, 'Diff');
    tTypes = rt(~isDiff);
    diffTypes = rt(isDiff);

    % Sort T entries numerically
    tVals = cellfun(@(s) str2double(regexp(s, '\d+', 'match', 'once')), tTypes);
    [~, tIdx] = sort(tVals);
    tTypes = tTypes(tIdx);

    % Parse the two T-values from each Diff entry
    vals = zeros(numel(diffTypes), 2);
    for i = 1:numel(diffTypes)
        tokens = regexp(diffTypes{i}, 'T(\d+)', 'tokens');
        vals(i,:) = sort([str2double(tokens{1}{1}), str2double(tokens{2}{1})]);
    end

    % Sort by span (sequential first, largest span last), then by smaller T value
    span = vals(:,2) - vals(:,1);
    [~, dIdx] = sortrows([span, vals(:,1)]);
    diffTypes = diffTypes(dIdx);

    % Combine: all T entries first, then all Diff entries
    ort = [tTypes(:); diffTypes(:)]';
end

function ort = get_orderedRowTypes(expinfo)
    rt = arrayfun(@(x) x.metrics.Type, expinfo, 'UniformOutput', false);
    rt = unique(vertcat(rt{:}), 'stable');
    % Extract only the Diff entries
    isDiff = startsWith(rt, 'Diff');
    diffTypes = rt(isDiff);
    
    % Parse the two T-values from each Diff entry
    vals = zeros(numel(diffTypes), 2);
    for i = 1:numel(diffTypes)
        tokens = regexp(diffTypes{i}, 'T(\d+)', 'tokens');
        vals(i,:) = sort([str2double(tokens{1}{1}), str2double(tokens{2}{1})]);
    end
    
    % Sort by span (sequential first, largest span last), then by smaller T value
    span = vals(:,2) - vals(:,1);
    [~, idx] = sortrows([span, vals(:,1)]);
    diffTypes = diffTypes(idx);
    vals = vals(idx,:);
    
    % Build the output: [smaller_T, larger_T, Diff] for each
    ort = cell(1, 3*numel(diffTypes));
    for i = 1:numel(diffTypes)
        ort(3*i-2 : 3*i) = {sprintf('T%d', vals(i,1)), sprintf('T%d', vals(i,2)), diffTypes{i}};
    end
end

function r = parse_cycleMetricsResults(results, halfCycle, fs)
    r.meanCen = results.eye.centroidMean(:);
    r.meanCenAdj = (r.meanCen + [0; halfCycle]) / fs;
    r.meanLag = results.centroidMeanDiff(:);
    r.meanSkew = results.eye.skewFromCentroidMean(:);

    r.medCen = results.eye.centroidMedian(:);
    r.medCenAdj = (r.medCen + [0; halfCycle]) / fs;
    r.medLag = results.centroidMedianDiff(:);
    r.medSkew = results.eye.skewFromCentroidMedian(:);

    r.bowleySkew = results.eye.skewQuantile(:);
    r.meanmedSkew = results.eye.skewMeanMedian(:);

    r.meanTxt = cellstr(compose("MEAN (dotted):\n    Centroid = %.1f ms\n    Lag = %.1f ms\n    Skew = %.3f\n    MeanMedSkew = %.3f", ...
        r.meanCen(:), r.meanLag(:), r.meanSkew(:), r.meanmedSkew(:)));
    r.medTxt = cellstr(compose("MEDIAN (dashed):\n    Centroid = %.1f ms\n    Lag = %.1f ms\n    Skew = %.3f\n    BowleySkew = %.3f", ...
        r.medCen(:), r.medLag(:), r.medSkew(:), r.bowleySkew(:)));

    % Find the two tallest peaks in the cycle-mean difference data
    r.ntPks = [results.eye.peak1Val(1), results.eye.peak2Val(1)];
    r.ntLocsMs = [results.eye.peak1TimeMs(1), results.eye.peak2TimeMs(1)];
    r.ntLocs = round(r.ntLocsMs) ./ fs;
    r.tnPks = [results.eye.peak1Val(2), results.eye.peak2Val(2)];
    r.tnLocsMs = [results.eye.peak1TimeMs(2), results.eye.peak2TimeMs(2)];
    r.tnLocs = round(r.tnLocsMs + round(halfCycle)) ./ fs;
end

function [expinfo,nexp] = get_expInfo(source_folder)
    % Recursively find all .xlsx files in source folder
    expinfo = dir(fullfile(source_folder, '**', '*_analysis-diffdata.xlsx'));
    
    % Remove unused fields
    expinfo = expinfo([expinfo.isdir]==0);
    nexp = length(expinfo);
    expinfo = rmfield(expinfo, {'date','bytes','datenum','isdir'});
    
    % Extract key-value metadata from filename and add as columns
    for i = 1:numel(expinfo)
        [~, fname] = fileparts(expinfo(i).name);
        tokens = regexp(fname, '(?<key>[^_-]+)-(?<val>[^_]+)', 'names');
        for j = 1:numel(tokens)
            expinfo(i).(tokens(j).key) = tokens(j).val;
        end
    end   
end
