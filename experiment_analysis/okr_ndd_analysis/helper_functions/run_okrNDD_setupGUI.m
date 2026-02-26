function [params,d] = run_okrNDD_setupGUI(params)
%RUN_OKRNDD_SETUPGUI GUI for user to set analysis parameters

% Set flag to check whether all three steps in setup was completed properly
params.setupAllCompleted = [0 0 0];

% Initialize other output variables
d = struct;
smr = [];
scaleCh1 = [];
scaleCh2 = [];

% Estimate the appropriate size for figure
figsize = calc_figureSize();


%% Initialize GUI main figure
fig = uifigure('Toolbar','none', ... 
    'MenuBar','none', ...
    'NumberTitle','off', ...
    'Name','OKR NDD Analysis App', ...
    'Position',figsize, ...
    'WindowState','maximized');

% Add confirmation dialog on close
fig.CloseRequestFcn = @(src,evt) on_Close(src);

% Define GUI parameters
fontsz1 = 12;


%% Main grid layout for entire figure
GridLayout_Main = uigridlayout(fig, ...
    'RowHeight',{80,50,'1x'}, ...  % Sets number of rows and their heights
    'ColumnWidth',{'1x'}, ...  % Sets number of columns and their widths
    'Padding',[20 20 20 20]);  % Sets padding from the edge of the figure


%% Sub-grid layout for setting experiment file, calibration scaling, and save folder
GridLayout_Paths = uigridlayout(GridLayout_Main, ...
    'RowHeight',{'1x', '1x', '1x'}, ...
    'ColumnWidth',{130, 80, 80, '1x', 200}, ...
    'Padding',[0 0 0 0]);
GridLayout_Paths.Layout.Row = 1;

% Component for loading experiment file
LabelLoadExp = uilabel(GridLayout_Paths, ...
    'Text','Experiment File:', ...
    'VerticalAlignment','center', ...
    'FontSize',fontsz1);
LabelLoadExp.Layout.Row = 1;
LabelLoadExp.Layout.Column = 1;
ButtonLoadExp = uibutton(GridLayout_Paths, ...
    'Text','Browse', ...
    'FontSize',fontsz1, ...
    'VerticalAlignment','center', ...
    'ButtonPushedFcn',@(src,evt) selectFile('*.mat;*.smr', 'exp_filepath'));
ButtonLoadExp.Layout.Row = 1;
ButtonLoadExp.Layout.Column = 2;
TextFieldExp = uieditfield(GridLayout_Paths, ...
    'InputType','text', ...
    'Placeholder',' Load Spike2 experiment recording (.smr or .mat) file to analyze.', ...
    'FontSize',fontsz1, ...
    'Editable','off');
TextFieldExp.Layout.Row = 1;
TextFieldExp.Layout.Column = [3 4];

% Component for loading calibration file
LabelLoadCal = uilabel(GridLayout_Paths, ...
    'Text','Calibration Scaling:', ...
    'FontSize',fontsz1);
LabelLoadCal.Layout.Row = 2;
LabelLoadCal.Layout.Column = 1;
ButtonLoadCal = uibutton(GridLayout_Paths, ...
    'Text','Browse', ...
    'FontSize',fontsz1, ...
    'VerticalAlignment','center', ...
    'ButtonPushedFcn',@(src,evt) selectFile('*cal*.mat', 'cal_filepath'));
ButtonLoadCal.Layout.Row = 2;
ButtonLoadCal.Layout.Column = 2;
ButtonCustomCal = uibutton(GridLayout_Paths, ...
    'Text','Custom', ...
    'FontSize',fontsz1, ...
    'VerticalAlignment','center', ...
    'ButtonPushedFcn',@(src,evt) setCustomCal());
ButtonCustomCal.Layout.Row = 2;
ButtonCustomCal.Layout.Column = 3;
TextFieldCal = uieditfield(GridLayout_Paths, ...
    'InputType','text', ...
    'Placeholder',' Load scaling factors from calibration .mat file OR provide custom values.', ...
    'FontSize',fontsz1, ...
    'Editable','off');
TextFieldCal.Layout.Row = 2;
TextFieldCal.Layout.Column = 4;

% Component for setting optional target save folder location
LabelSaveDir = uilabel(GridLayout_Paths, ...
    'Text','Save Analysis to:', ...
    'VerticalAlignment','center', ...
    'FontSize',fontsz1);
LabelSaveDir.Layout.Row = 3;
LabelSaveDir.Layout.Column = 1;
ButtonSaveDir = uibutton(GridLayout_Paths, ...
    'Text','Browse', ...
    'FontSize',fontsz1, ...
    'VerticalAlignment','center', ...
    'ButtonPushedFcn',@(src,evt) selectFolder('save_folderpath'));
ButtonSaveDir.Layout.Row = 3;
ButtonSaveDir.Layout.Column = 2;
saveDirText = [' (Optional) Set folder to save analysis results in.', ...
    ' If left empty, analysis will automatically be saved in the experiment folder.'];
TextFieldSave = uieditfield(GridLayout_Paths, ...
    'InputType','text', ...
    'Placeholder',saveDirText, ...
    'FontSize',fontsz1, ...
    'Editable','off');
TextFieldSave.Layout.Row = 3;
TextFieldSave.Layout.Column = [3 4];

% Create next button
ButtonNext = uibutton(GridLayout_Paths, ...
    'Text','NEXT', ...
    'FontWeight','bold', ...
    'FontSize',18, ...
    'Enable','off', ...
    'ButtonPushedFcn',@(src,evt) nextButtonPushed());
ButtonNext.Layout.Row = [1 3];
ButtonNext.Layout.Column = 5;


%% Sub-grid panel and layout for setting saccade threshold interactively
Panel_Params = uipanel(GridLayout_Main);
Panel_Params.Layout.Row = 2;
columnwidths = horzcat({'0.25x',105,75}, ...
    repmat({'0.25x',105,75}, 1, 6), ...
    {'1x'});
GridLayout_Params = uigridlayout(Panel_Params, ...
    'RowHeight',{'1x'}, ...
    'ColumnWidth',columnwidths, ...
    'ColumnSpacing',0); 

% Left block number label and dropdown components
LabelLeftBlock = uilabel(GridLayout_Params, ...
    'Text','Left Block: ', ...
    'FontSize',fontsz1, ...
    'HorizontalAlignment','right');
LabelLeftBlock.Layout.Column = 2;
DropDownLeftBlock = uidropdown(GridLayout_Params, ...
    'Items',{'None'}, ...
    'ValueIndex',1, ...
    'FontSize',fontsz1, ...
    'ValueChangedFcn',@(src,evt) dropdownChanged(evt, 2, [1,2,3]));

LabelCohort = uilabel(GridLayout_Params, ...
    'Text','Cohort: ', ...
    'FontSize',fontsz1, ...
    'HorizontalAlignment','right');
LabelCohort.Layout.Column = 5;
TextFieldCohort = uieditfield(GridLayout_Params, ...
    'InputType','alphanumerics', ...
    'FontSize',fontsz1, ...
    'Value','None', ...
    'HorizontalAlignment','center');

% Text components for subject ID and cohort
LabelSubId = uilabel(GridLayout_Params, ...
    'Text','Mouse ID: ', ...
    'FontSize',fontsz1, ...
    'HorizontalAlignment','right');
LabelSubId.Layout.Column = 8;
TextFieldSubID = uieditfield(GridLayout_Params, ...
    'InputType','alphanumerics', ...
    'FontSize',fontsz1, ...
    'Value','None', ...
    'HorizontalAlignment','center');

% Subject and task condition label and dropdown components
LabelSubCond = uilabel(GridLayout_Params, ...
    'Text','Mouse Condition: ', ...
    'FontSize',fontsz1, ...
    'HorizontalAlignment','right');
LabelSubCond.Layout.Column = 11;
DropDownSubCond = uidropdown(GridLayout_Params, ...
    'Items',params.subcond_options, ...
    'ValueIndex',1, ...
    'FontSize',fontsz1);

LabelTaskCond = uilabel(GridLayout_Params, ...
    'Text','Task Condition: ', ...
    'FontSize',fontsz1, ...
    'HorizontalAlignment','right');
LabelTaskCond.Layout.Column = 14;
DropDownTaskCond = uidropdown(GridLayout_Params, ...
    'Items',params.taskcond_options, ...
    'ValueIndex',1, ...
    'FontSize',fontsz1);

% Saccade threshold label and numeric field components
LabelThresh = uilabel(GridLayout_Params, ...
    'Text','Saccade Thresh: ', ...
    'FontSize',fontsz1, ...
    'HorizontalAlignment','right');
LabelThresh.Layout.Column = 17;
NumericFieldThresh = uieditfield(GridLayout_Params, ...
    'numeric', ...
    'FontSize',fontsz1, ...
    'Value',params.saccadeThresh, ...
    'HorizontalAlignment','center', ...
    'ValueChangedFcn',@(src,evt) threshChanged(evt));

% Right block number label and dropdown components
LabelRightBlock = uilabel(GridLayout_Params, ...
    'Text','Right Block: ', ...
    'FontSize',fontsz1, ...
    'HorizontalAlignment','right');
LabelRightBlock.Layout.Column = 20;
DropDownRightBlock = uidropdown(GridLayout_Params, ...
    'Items',{'None'}, ...
    'ValueIndex',1, ...
    'FontSize',fontsz1, ...
    'ValueChangedFcn',@(src,evt) dropdownChanged(evt, 3, [4,5,6]));


%% Sub-grid layout for displaying saccade threshold results subplots
GridLayout_Plots = uigridlayout(GridLayout_Main, ...
    'RowHeight',{'0.4x', '1x'}, ...
    'ColumnWidth',{'1x', '1x'});
GridLayout_Plots.Layout.Row = 3;

% Define line colors for plots
c1 = 'k';
c2 = 'b';
c3 = 'r';
ylims = [-100, 100];

% Initialize uiaxes for nGoodCycles bar subplot
axs(1) = uiaxes(GridLayout_Plots);
axesCustomToolbarButtons(axs(1));
axs(1).Layout.Row = 1;
axs(1).Layout.Column = [1 2];
title(axs(1), 'Number of Good Cycles per Block');
hold(axs(1), 'on');
axbar = bar(axs(1), NaN, NaN);
axtext = text(axs(1), NaN, NaN, '', 'HorizontalAlignment','center', 'VerticalAlignment','bottom');
hold(axs(1), 'off');
ylim(axs(1), [0 1]);
xlabel(axs(1), 'Block Number', 'FontSize',fontsz1);
ylabel(axs(1), 'Good Cycle Fraction', 'FontSize',fontsz1);
box(axs(1), 'on');
grid(axs(1), 'on');
axs(1).GridAlpha = 0.06;

% Initialize uiaxes for left block subplot
axs(2) = uiaxes(GridLayout_Plots);
axesCustomToolbarButtons(axs(2));
axs(2).Layout.Row = 2;
axs(2).Layout.Column = 1;
title(axs(2), 'Block NaN (NaN of NaN): ampSEM = NaN | phaseSEM = NaN  | varRes = NaN');
hold(axs(2), 'on');
yline(axs(2), 0, 'Color',[0 0 0]+0.7, 'LineWidth',0.1, ...
    'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
axlines(1) = plot(axs(2), NaN, NaN, ...
    'Color',c1, 'LineWidth',1, 'LineStyle','--', ...
    'DisplayName',' Stimulus Velocity', ...
    'HitTest','off', 'PickableParts','none');
axlines(2) = plot(axs(2), NaN, NaN, ...
    'Color',c3, 'LineWidth',1, ...
    'DisplayName',' Saccades', ...
    'HitTest','off', 'PickableParts','none');
axlines(3) = plot(axs(2), NaN, NaN, ...
    'Color',c2, 'LineWidth',1, ...
    'DisplayName',' Desaccaded Eye Velocity', ...
    'HitTest','off', 'PickableParts','none'); 
axlines(4) = plot(axs(2), NaN, NaN, ...
    'Color',c1, 'LineWidth',2, ...
    'DisplayName',' Fit of Raw Eye Velocity', ...
    'HitTest','off', 'PickableParts','none');
hold(axs(2), 'off');
xlabel(axs(2), 'Time (s)', 'FontSize',fontsz1);
ylabel(axs(2), 'Velocity (deg/s)', 'FontSize',fontsz1);
ylim(axs(2), ylims);
box(axs(2), 'on');
grid(axs(2), 'on');
axs(2).GridAlpha = 0.06;
legend(axs(2), 'Location','northwest', 'FontSize',fontsz1-2);

% Initialize uiaxes for right block subplot
axs(3) = uiaxes(GridLayout_Plots);
axesCustomToolbarButtons(axs(3));
axs(3).Layout.Row = 2;
axs(3).Layout.Column = 2;
title(axs(3), 'Block NaN (NaN of NaN): ampSEM = NaN | phaseSEM = NaN  | varRes = NaN');
hold(axs(3), 'on');
yline(axs(3), 0, 'Color',[0 0 0]+0.7, 'LineWidth',0.1, ...
    'HandleVisibility','off', 'HitTest','off', 'PickableParts','none');
axlines(5) = plot(axs(3), NaN, NaN, ...
    'Color',c1, 'LineWidth',1, 'LineStyle','--', ...
    'DisplayName',' Stimulus Velocity', ...
    'HitTest','off', 'PickableParts','none');
axlines(6) = plot(axs(3), NaN, NaN, ...
    'Color',c3, 'LineWidth',1, ...
    'DisplayName',' Saccades', ...
    'HitTest','off', 'PickableParts','none');
axlines(7) = plot(axs(3), NaN, NaN, ...
    'Color',c2, 'LineWidth',1, ...
    'DisplayName',' Desaccaded Eye Velocity', ...
    'HitTest','off', 'PickableParts','none'); 
axlines(8) = plot(axs(3), NaN, NaN, ...
    'Color',c1, 'LineWidth',2, ...
    'DisplayName',' Fit of Raw Eye Velocity', ...
    'HitTest','off', 'PickableParts','none');
hold(axs(3), 'off');
xlabel(axs(3), 'Time (s)', 'FontSize',fontsz1);
ylim(axs(3), ylims);
box(axs(3), 'on');
grid(axs(3), 'on');
axs(3).GridAlpha = 0.06;
legend(axs(3), 'Location','northwest', 'FontSize',fontsz1-2);
linkaxes([axs(2), axs(3)], 'xy');


uiwait(fig);


    %% Helper functions
    function titleTxt = get_blockTitleText(bidx, ngc, ntc)
        cvd = d.all_cvd{bidx};
        titleTxt = [sprintf('Block %d (%d of %d):', bidx, ngc, ntc), ...
            sprintf(' ampSEM = %.4f | phaseSEM = %.4f  | varRes = %.4f', cvd.amplitudeSEM, cvd.phaseSEM_deg, cvd.varResidualMean)];
    end

    function updateBlockPlot(blockId, axIdx, lineIds)
        stimvel_ii = d.all_stimvel{blockId};
        hevel_ii = d.all_hevel{blockId};
        hevel_des_ii = d.all_hevel_des{blockId};
        hevel_rawfit_ii = d.all_hevel_rawfit{blockId};
        btimes = (0:length(hevel_ii)-1) / params.fs;
        btxt = get_blockTitleText(blockId, ...
            d.all_nGoodCycles(blockId), ...
            d.all_nTotalCycles(blockId));
        title(axs(axIdx), btxt);
        axlines(lineIds(1)).XData = btimes;
        axlines(lineIds(1)).YData = stimvel_ii;
        axlines(lineIds(2)).XData = btimes;
        axlines(lineIds(2)).YData = hevel_ii;
        axlines(lineIds(3)).XData = btimes;
        axlines(lineIds(3)).YData = hevel_des_ii;
        axlines(lineIds(4)).XData = btimes;
        axlines(lineIds(4)).YData = hevel_rawfit_ii;
        xlim(axs(axIdx), [0 btimes(end)]);
        drawnow limitrate;
    end

    function updateBarPlot()
        blockNumbers = 1:d.nBlocks;
        % Build the full CData matrix before assigning
        cdata = zeros(d.nBlocks, 3);
        cdata(d.all_nGoodCycles > 15, :)  = repmat([0 0.7 0], sum(d.all_nGoodCycles > 15), 1);
        cdata(d.all_nGoodCycles <= 15, :) = repmat([1 0.5 0], sum(d.all_nGoodCycles <= 15), 1);
        cdata(d.all_nGoodCycles < 11, :)  = repmat([0.7 0 0], sum(d.all_nGoodCycles < 11), 1);
        % Update bar data and colors together
        set(axbar, 'XData',blockNumbers, 'YData',d.all_goodCyclesFrac(:), ...
            'CData',cdata, 'FaceColor','flat');
        blabels = string(d.all_nGoodCycles);
        delete(axtext);
        axtext = text(axs(1), axbar.XEndPoints, axbar.YEndPoints, blabels, ...
            'FontSize',fontsz1-2, ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom');
        % Update bar plot xlim
        xlim(axs(1), [0.4 d.nBlocks+0.6]);
        xticks(axs(1), blockNumbers);
        xticklabels(axs(1), cellstr(string(blockNumbers)));
        drawnow limitrate;
    end

    function reprocessSaccadeThreshold()
        nBlocks = d.nBlocks;
        thresh = NumericFieldThresh.Value;
    
        all_nGoodCycles = nan(nBlocks, 1);
        all_nTotalCycles = nan(nBlocks, 1);
        all_goodCyclesFrac = nan(nBlocks, 1);
        all_hevel_des = cell(nBlocks, 1);
        all_hevel_goodcycles = cell(nBlocks, 1);
        all_cvd = cell(nBlocks, 1);
    
        for ii = 1:nBlocks
            hevel_ii       = d.all_hevel{ii};
            hevel_raw_fit  = d.all_hevel_rawfit{ii};
            startpt        = d.all_startpt(ii);
            hevel_cyclemat = d.all_hevel_cyclemat{ii};
    
            % Recompute saccade mask with new threshold
            [~,saccMask,hevel_des_ii] = run_saccadeDetection(hevel_ii, ...
                hevel_raw_fit, ...
                thresh, ...
                d.saccadePad, ...
                params.minGoodChunk_len, ...
                params.saccadeMethod);
    
            % Segment and find bad cycles
            sacc_mat = segmentIntoCycles(double(saccMask), startpt, d.cycleLength);
            badCycles = any(sacc_mat, 2);
    
            all_nGoodCycles(ii)      = sum(~badCycles);
            [all_nTotalCycles(ii),~] = size(sacc_mat);
            all_goodCyclesFrac(ii)   = all_nGoodCycles(ii) / all_nTotalCycles(ii);
    
            if all_nGoodCycles(ii) > 0
                hevel_good_cyclemat = hevel_cyclemat(~badCycles, :);
            else
                hevel_good_cyclemat = zeros(size(hevel_cyclemat));
            end
    
            all_hevel_des{ii} = hevel_des_ii;
            all_hevel_goodcycles{ii} = hevel_good_cyclemat;
            all_cvd{ii} = cycleVarianceDecomposition(hevel_good_cyclemat);
        end
    
        % Update data struct (only the fields that changed)
        d.all_nGoodCycles = all_nGoodCycles;
        d.all_nTotalCycles = all_nTotalCycles;
        d.all_goodCyclesFrac = all_goodCyclesFrac;
        d.all_hevel_des = all_hevel_des;
        d.all_hevel_goodcycles = all_hevel_goodcycles;
        d.all_cvd = all_cvd;
    
        % Refresh plots
        updateBarPlot();
        updateBlockPlot(DropDownLeftBlock.ValueIndex, 2, [1,2,3,4]);
        updateBlockPlot(DropDownRightBlock.ValueIndex, 3, [5,6,7,8]);
    end

    function processData()
        % Display progress dialog
        pdlg = uiprogressdlg(fig, 'Title','Importing data, please wait...');
        
        % Validate the value of nBlocks
        nBlocks = d.nBlocks;
        if isempty(nBlocks) || ~isnumeric(nBlocks) || (nBlocks < 1)
            error('Invalid value of nBlocks: %s! Abortin script...', string(nBlocks));
        end
        
        % Validate the blockIds array
        blockIds = d.blockIds;
        if isempty(blockIds) || ~isequal(size(blockIds,1), 2) || ~isequal(size(blockIds,2), nBlocks)
            error('The array "blockIds" is invalid! Aborting script...');
        end

        all_nGoodCycles = nan(nBlocks, 1); 
        all_nTotalCycles = nan(nBlocks, 1); 
        all_goodCyclesFrac = nan(nBlocks, 1); 
        all_startpt = nan(nBlocks, 1);
        all_stimvel = cell(nBlocks, 1); 
        all_hevel_cyclemat = cell(nBlocks, 1); 
        all_hevel = cell(nBlocks, 1); 
        all_hevel_des = cell(nBlocks, 1); 
        all_hevel_rawfit = cell(nBlocks, 1); 
        all_hevel_goodcycles = cell(nBlocks, 1);
        all_cvd = cell(nBlocks, 1);

        for ii = 1:nBlocks
            % Extract corresponding block segments for reference stimulus (chair/drum) and eye data
            chairvel_raw_ii = d.chairvel_raw(blockIds(1,ii):blockIds(2,ii));
            drumvel_raw_ii  = d.drumvel_raw(blockIds(1,ii):blockIds(2,ii));
            hevel_raw_ii = d.hevel_raw(blockIds(1,ii):blockIds(2,ii));
            hepos_raw_ii = d.hepos_raw(blockIds(1,ii):blockIds(2,ii));
            % Get filtered eye position and eye velocity traces
            hepos_ii = butterworthfilter(hepos_raw_ii, params.lowpassCutoff, params.fs, 'n',4);
            
            % % TODO: Impliment FIF?
            % opts = Settings_FIF_v3('delta',params.fifSettings.delta, ...
            %     'Xi',params.fifSettings.Xi, ...
            %     'NIMFs',params.fifSettings.NIMFs, ...
            %     'alpha',params.fifSettings.alpha, ...
            %     'verbose',0);
            % [IMFs,~] = FIF_v2_12(hepos_raw_ii(:)', opts);
            % hepos_ii = sum(IMFs(params.fifSettings.minIMF,:), 1);

            hevel_ii = movingslope(hepos_ii, params.filterWindow) * params.fs;
            % Perform initial fit of the raw eye velocity data
            keep = abs(hevel_raw_ii) < 5 * std(abs(hevel_raw_ii)) + mean(abs(hevel_raw_ii));
            [fit0,~] = calc_sineFit(params.exp_stimfreq, params.fs, keep, hevel_raw_ii);
            hevel_raw_fit = fit0.eyevel_fit;
            % Apply saccade detection to the filtered eye velocity trace to get the necessary saccade masks
            [saccDilMask,saccMask,hevel_des_ii] = run_saccadeDetection(hevel_ii, ...
                hevel_raw_fit, ...
                NumericFieldThresh.Value, ...
                d.saccadePad, ...
                params.minGoodChunk_len, ...
                params.saccadeMethod);
            % Calculate sinusoidal fits of the desaccaded eye, chair, and drum velocity block data
            [bfit,~] = calc_sineFit(params.exp_stimfreq, params.fs, ~saccDilMask, hevel_ii, chairvel_raw_ii, drumvel_raw_ii);
            % Calculate the start point of the first positive stimulus cycle
            startpt = max(1, round(mod(-bfit.ref_angle,360) / 360 * params.fs/params.exp_stimfreq));
            % Get cycles
            hevel_cyclemat = segmentIntoCycles(hevel_ii, startpt, d.cycleLength); 
            sacc_mat       = segmentIntoCycles(double(saccMask), startpt, d.cycleLength);
            % sacc_mat is used to find the "bad" cycles that contain any NaN values
            badCycles = any(sacc_mat, 2);
            % This is used to calculate how many "good" cycles there are
            all_nGoodCycles(ii)      = sum(~badCycles);
            [all_nTotalCycles(ii),~] = size(sacc_mat);
            all_goodCyclesFrac(ii)   = all_nGoodCycles(ii) / all_nTotalCycles(ii);
            % A separate cycle mean is calculated using only "good" cycles 
            if all_nGoodCycles(ii) > 0
                hevel_good_cyclemat = hevel_cyclemat(~badCycles, :);
            else
                hevel_good_cyclemat = zeros(size(hevel_cyclemat));
            end
            all_startpt(ii) = startpt;
            if contains(string(d.blockTypes(ii)), 'OKR')
                all_stimvel{ii} = drumvel_raw_ii;
            else
                all_stimvel{ii} = chairvel_raw_ii;
            end
            all_hevel_cyclemat{ii} = hevel_cyclemat;
            all_hevel{ii} = hevel_ii;
            all_hevel_des{ii} = hevel_des_ii;
            all_hevel_rawfit{ii} = hevel_raw_fit;
            all_hevel_goodcycles{ii} = hevel_good_cyclemat;
            all_cvd{ii} = cycleVarianceDecomposition(hevel_good_cyclemat);
        end
        
        % Update data struct
        d.all_goodCyclesFrac = all_goodCyclesFrac;
        d.all_nGoodCycles = all_nGoodCycles;
        d.all_nTotalCycles = all_nTotalCycles;
        d.all_startpt = all_startpt;
        d.all_stimvel = all_stimvel;
        d.all_hevel_cyclemat = all_hevel_cyclemat;
        d.all_hevel = all_hevel;
        d.all_hevel_des = all_hevel_des;
        d.all_hevel_rawfit = all_hevel_rawfit;
        d.all_hevel_goodcycles = all_hevel_goodcycles;
        d.all_cvd = all_cvd;

        % Update the bar plot
        updateBarPlot();

        % Update selected left/right block dropdowns
        DropDownLeftBlock.ValueIndex = params.timepoint_ids(1);
        DropDownRightBlock.ValueIndex = params.timepoint_ids(end);

        % Update the respective plots
        updateBlockPlot(params.timepoint_ids(1), 2, [1,2,3,4]);
        updateBlockPlot(params.timepoint_ids(end), 3, [5,6,7,8]);
        
        % Close progress dialog
        close(pdlg);
    end

    function updatePanelParams()
        % Get and update block numbers
        if isfield(d, 'blockTypes')
            blockNumbers = cellstr(string(1:length(d.blockTypes)));
            DropDownLeftBlock.Items = blockNumbers;
            DropDownRightBlock.Items = blockNumbers;
        end
        % Attempt to extract required metadata from filename or prompt user
        % Subject ID
        if isfield(params, 'exp_name')
            params.exp_subid = regexp(params.exp_name, 'WT[-_]?\d+(?=[_\- ])', 'match', 'once');
            params.exp_subid = erase(params.exp_subid, {'-', '_', ' '});
        end
        if ~isfield(params, 'exp_subid') || isempty(params.exp_subid)
            TextFieldSubID.Value = 'None';
        else
            TextFieldSubID.Value = params.exp_subid;
        end
        % Subject condition
        if isfield(params, 'exp_filename')
            subcond_matches = find(cellfun(@(x) ~isempty(regexpi(params.exp_filename, x, 'once')), params.subcond_options), 1);
            if isequal(length(subcond_matches), 1)
                params.exp_subcond = params.subcond_options{subcond_matches};
                DropDownSubCond.ValueIndex = subcond_matches;
            else
                params.exp_subcond = params.subcond_options{1};
                DropDownSubCond.ValueIndex = 1;
            end
        end
    end

    function getExpData(smr, scaleCh1, scaleCh2)
        % Display progress dialog
        pdlg = uiprogressdlg(fig, 'Title','Importing data, please wait...');
        % Retrieve main data needed for analysis
        fprintf('\nImporting data, please wait:\n');
        % Extract the start/end times of each block and their types
        disp('    Identifying the start and end times of experiment blocks...');
        kbidx = find(strcmp(smr.channeltitles,'Keyboard'));  % Index of Keyboard channel in smr
        if isempty(kbidx)
            error('Spike2 file missing required "Keyboard" channel! Aborting script...\n'); 
        end
        kbmd = cell2mat(smr.channels{kbidx}.markerdata);     % Keyboard marker channel data
        kbmt = smr.channels{kbidx}.markertimes;              % Keyboard marker timestamps
        bss  = blockSchemaSearch(kbmd, [], []);              % Identify the start/end of blocks
        d.nBlocks    = bss.nall;                   % Number of blocks in the experiment
        d.blockTypes = categorical(bss.types, ...  % Assign detected type of each block
            [1 2 3 4 0], {'VORD', 'OKR', 'VORD', 'VORD', 'ERROR'});
        d.bssTimes = kbmt(bss.all);              % Start/end times in seconds of each block
        % The raw hepos data is finally loaded and converted to double precision
        disp('    Loading raw eye, chair, and drum data...');
        eye1idx = find(strcmp(smr.channeltitles,'hepos1'));
        eye2idx = find(strcmp(smr.channeltitles,'hepos2'));
        if isempty(eye1idx) || isempty(eye2idx) 
            error('Spike2 file missing required "hepos1"/"hepos2" channel! Aborting script...\n'); 
        end
        hepos1 = smr.channels{eye1idx,1}.data;
        hepos2 = smr.channels{eye2idx,1}.data;
        hepos_raw = double(scaleCh1*hepos1 + scaleCh2*hepos2);
        % Remove transients from the raw eye position data
        d.hepos_raw = removeTransients(hepos_raw, params.transientThresh);  % BRIAN
        % Load the raw chair and drum velocity data (applying transients removal too)
        % Note: These two lines of code aren't explicitly performed in the original pipeline
        chairidx = find(strcmp(smr.channeltitles,'hhvel'));
        drumidx = find(strcmp(smr.channeltitles,'htvel'));
        if isempty(chairidx) || isempty(drumidx)
            error('Spike2 file missing required "hhvel"/"htvel" channel! Aborting script...\n'); 
        end
        d.chairvel_raw = removeTransients(double(smr.channels{chairidx,1}.data), params.transientThresh);
        d.drumvel_raw  = removeTransients(double(smr.channels{drumidx,1}.data), params.transientThresh);
        % Extract the data samplerate from the experiment file
        params.fs = smr.channels{drumidx,1}.samplerate;
        % The raw hevel data is derived via SG filter
        % Note: The original code used a hardcoded 10-ms window!
        d.hevel_raw = movingslope(d.hepos_raw, params.filterWindow) * params.fs;
        % Convert blockTimes into respective start/end indices of each block
        d.blockIds = min(max(1, round(d.bssTimes * params.fs)), length(d.drumvel_raw));
        % Convert saccade detection padding to ms
        d.saccadePad = round(params.saccadeLRPad * params.fs);
        % Initialize length and time vectors for cycles (always the same)
        d.cycleLength = round(params.fs / params.exp_stimfreq);
        d.cycleTimes  = (1:d.cycleLength) / params.fs;
        % Update parameter panel
        updatePanelParams();
        % Close progress dialog
        close(pdlg);
        % Process data with updated params
        processData();
    end

    function loadExpFile(filepath)
        % If calibration was already loaded, ask user to confirm it still applies
        if params.setupAllCompleted(2)
            selection = uiconfirm(fig, ...
                sprintf('Keep current calibration (scaleCh1=%.4f, scaleCh2=%.4f) for the new file?', scaleCh1, scaleCh2), ...
                'Calibration Check', ...
                'Options', {'Keep', 'Clear'}, ...
                'DefaultOption', 'Clear');
            if strcmp(selection, 'Clear')
                % Reset GUI variables
                scaleCh1 = [];
                scaleCh2 = [];
                params.setupAllCompleted(2) = 0;
                params.cal_filepath = '';
                params.customScaleChs = [];
                TextFieldCal.Value = '';
                ButtonNext.Enable = 'off';
                % Reset plots to initial empty state
                for k = 1:6
                    axlines(k).XData = NaN;
                    axlines(k).YData = NaN;
                end
                title(axs(2), 'Block NaN (NaN of NaN): ampSEM = NaN | phaseSEM = NaN  | varRes = NaN');
                title(axs(3), 'Block NaN (NaN of NaN): ampSEM = NaN | phaseSEM = NaN  | varRes = NaN');
                % Reset bar plot
                set(axbar, 'XData', NaN, 'YData', NaN, 'CData',[0 0 0]);
                delete(axtext);
                axtext = text(axs(1), NaN, NaN, '', 'HorizontalAlignment','center', 'VerticalAlignment','bottom');
                % Reset dropdowns
                DropDownLeftBlock.Items = {'None'};
                DropDownRightBlock.Items = {'None'};
                drawnow limitrate;
            end
        end

        % Display progress dialog
        pdlg = uiprogressdlg(fig, 'Title','Loading file, please wait...');
        % Try to load the experiment file and update text field if successful
        try
            [exp_folderpath,exp_name,exp_ext] = fileparts(filepath);
            exp_filename = [exp_name, exp_ext];
            % Load the experiment file 
            if strcmp(exp_ext, '.mat')
                temp = load(filepath).smr;
            elseif strcmp(exp_ext, '.smr')
                temp = importSpike2File(filepath);
            else
                TextFieldExp.Value = '';
                params.setupAllCompleted(1) = 0;
                ButtonNext.Enable = 'off';
                warning('Invalid Experiment file: %s', exp_filename);
                return;
            end
            % Verify imported data is valid
            if ~isfield(temp, 'filename')
                TextFieldExp.Value = '';
                params.setupAllCompleted(1) = 0;
                ButtonNext.Enable = 'off';
                warning('Invalid Experiment file: %s', exp_filename);
                return;
            end
            smr = temp;
            % Set text field and smr property
            TextFieldExp.Value = filepath;
            TextFieldExp.Tooltip = filepath;
            % Update params file with experiment file info
            params.exp_filepath   = filepath;
            params.exp_folderpath = exp_folderpath;
            params.exp_filename   = exp_filename;
            params.exp_name       = exp_name;
            params.exp_ext        = exp_ext;
            % Set first flag as completed
            params.setupAllCompleted(1) = 1;
            fprintf('Loaded experiment file: %s\n', exp_filename);
            % Close progress dialog
            close(pdlg);
            if all(params.setupAllCompleted(1:2))
                ButtonNext.Enable = 'on';
                getExpData(smr, scaleCh1, scaleCh2);
            else
                ButtonNext.Enable = 'off';
            end
        catch ME
            % Close progress dialog
            close(pdlg);
            TextFieldExp.Value = '';
            params.setupAllCompleted(1) = 0;
            ButtonNext.Enable = 'off';
            warning(ME.identifier, 'Error loading file: %s', ME.message);  
        end
    end

    function loadCalFile(filepath)
        try
            % Try to load the calibration file and update text field if successful
            % Verify calibration file contains valid scaling factors
            cal = load(filepath, 'scaleCh1', 'scaleCh2');
            % Verify calibration file has both scaling factors
            status1 = isfield(cal,'scaleCh1') && isfield(cal,'scaleCh2');
            % Verify only one scaling factor is nonzero
            try
                status2 = isequal(nnz([cal.scaleCh1,cal.scaleCh2]), 1);
            catch
                status2 = 0;
            end
            if ~status1 || ~status2
                TextFieldCal.Value = '';
                params.setupAllCompleted(2) = 0;
                ButtonNext.Enable = 'off';
                warning('Calibration .mat file has invalid scaling factors: %s', filepath);
                return;
            end
            scaleCh1 = cal.scaleCh1;
            scaleCh2 = cal.scaleCh2;
            % Set text field and scaleCh property
            TextFieldCal.Value = sprintf('scaleCh1 = %.4f | scaleCh2 = %.4f', cal.scaleCh1, cal.scaleCh2);
            % Update params file with calibration file info
            params.cal_filepath = filepath;
            % Set second flag as completed
            params.setupAllCompleted(2) = 1;
            fprintf('Loaded calibration scaling factors from file: %s\n', TextFieldCal.Value);
            if all(params.setupAllCompleted(1:2))
                ButtonNext.Enable = 'on';
                getExpData(smr, scaleCh1, scaleCh2);
            else
                ButtonNext.Enable = 'off';
            end
        catch ME
            TextFieldCal.Value = '';
            params.setupAllCompleted(2) = 0;
            ButtonNext.Enable = 'off';
            warning(ME.identifier, 'Error loading file: %s', ME.message);  
        end
    end

    function setCustomCal()
        % Give dialog to provide custom channel scaling values
        prompt = {'scaleCh1:', 'scaleCh2:'};
        dlgtitle = 'Custom Calibration Scaling Factors';
        fieldsize = [1 45; 1 45];
        defaults = {'', ''};
        errorTxt = ['Invalid value(s) was/were provided!', ...
            ' Both values must be provided.', ...
            ' Only one value can be nonzero.'];
        answer = inputdlg(prompt, dlgtitle, fieldsize, defaults);
        figure(fig);
        if isempty(answer)
            return % User cancelled
        end
        temp1 = str2double(answer{1});
        temp2 = str2double(answer{2});
        % Verify the values are valid
        if any(isnan([temp1,temp2])) || ~isequal(nnz([temp1,temp2]), 1)
            TextFieldCal.Value = '';
            params.setupAllCompleted(2) = 0;
            ButtonNext.Enable = 'off';
            uialert(fig, errorTxt, 'Invalid Input');
            return
        end
        scaleCh1 = temp1;
        scaleCh2 = temp2;
        % Set text field and scaleCh property
        TextFieldCal.Value = sprintf('scaleCh1 = %.4f | scaleCh2 = %.4f', temp1, temp2);
        % Update params file with custom scaling values
        params.customScaleChs = [temp1, temp2];
        % Set second flag as completed
        params.setupAllCompleted(2) = 1;
        fprintf('Custom calibration scaling factors provided: %s\n', TextFieldCal.Value);
        if all(params.setupAllCompleted(1:2))
            ButtonNext.Enable = 'on';
            getExpData(smr, scaleCh1, scaleCh2);
        else
            ButtonNext.Enable = 'off';
        end
    end

    function nout = formatName(sid, econd, etask, scohort)
        % Define naming format for all saved files and their respective folders.
        % This will make it easier for our group analysis scripts to recursively 
        % find all the relevant analysis data.
        % 
        % ename - Experiment file name (without the extension)
        % econd - Experiment subject condition (e.g., DR, NR, etc.)
        % etask - Experiment task condition (e.g. std, dim, dmr)
        % scohort - Subject cohort (optional, default [] or '')
        if ~exist('scohort','var') || isempty(scohort)
            scohort = "NONE";
        end
        nout = sprintf('cond-%s_task-%s_sub-%s_cohort-%s_analysis-diffdata', econd, etask, sid, scohort);
    end

    function nextButtonPushed()
        if isequal(params.setupAllCompleted(3), 0)
            params.save_folderpath = params.exp_folderpath;
            % Mark setup as fully completed and release the GUI
            params.setupAllCompleted(3) = 1;
        end
        % Save the currently set values on parameter panel
        params.exp_subcohort = TextFieldCohort.Value;
        params.exp_subid = TextFieldSubID.Value;
        params.exp_subcond = params.subcond_options{DropDownSubCond.ValueIndex};
        params.exp_taskcond = params.taskcond_options{DropDownTaskCond.ValueIndex};
        params.saccadeThresh = NumericFieldThresh.Value;
        % Define corresponding save folderpath
        save_foldername = formatName(params.exp_subid, params.exp_subcond, params.exp_taskcond, params.exp_subcohort);
        params.save_folderpath = fullfile(params.save_folderpath, save_foldername);
        % Prompt overwrite if the folderpath already exists
        if isfolder(params.save_folderpath)  % If so, ask to overwrite
            answer = questdlg('Save folder already exists, do you wish to overwrite previous results?', ...
                'Option to overwrite', 'Yes', 'Cancel', 'Yes');
            if strcmp(answer, 'Cancel')
                warning('Okay, analysis will not be saved!');
                params.save_folderpath      = 0;
                params.save_filepath        = 0;
                params.setupAllCompleted(3) = 0;
                return  % Do nothing
            end
        else  % If not, create the directory
            [~,~] = mkdir(params.save_folderpath);
        end
        % Define corresponding save filepath
        save_filename = [save_foldername, '.mat'];
        params.save_filepath = fullfile(params.save_folderpath, save_filename);
        fprintf('    Analysis results will be saved to: %s\n', save_filename);
        uiresume(fig);
        delete(fig);
    end

    function threshChanged(evt)
        NewValue = evt.Value;
        % Only continue if a valid value was given
        if isempty(NewValue) || isnan(NewValue) || (NewValue <= 0)
            warning('Provided saccade threshold is invalid and will be ignored: %s', string(NewValue));
            return
        end
        % Only reprocess if data has been loaded
        if all(params.setupAllCompleted(1:2))
            reprocessSaccadeThreshold();
        end
    end

    function dropdownChanged(evt, axidx, lineids)
        NewValue = evt.ValueIndex;
        updateBlockPlot(NewValue, axidx, lineids);
    end

    function selectFile(fileExt, paramsField)
        % Set main figure invisible to prevent file dialog being hidden
        fig.Visible = 'off';
        if isfield(params, 'exp_folderpath') && ~isempty(params.exp_folderpath)
            openLoc = params.exp_folderpath;
        else
            openLoc = pwd;
        end
        if strcmp(paramsField, 'exp_filepath')
            disp('    SELECT EXPERIMENT FILE (.SMR/.MAT) TO ANALYZE...');
        elseif strcmp(paramsField, 'cal_filepath')
            disp('    SELECT CALIBRATION FILE (.MAT) TO EXTRACT CHANNEL SCALING FACTORS...');
        end
        [filename,folderpath] = uigetfile(fileExt, ['Select ', fileExt], openLoc);
        fig.Visible = 'on';
        if isequal(filename, 0), return; end  % User canceled selection
        filepath = fullfile(folderpath, filename);
        params.(paramsField) = filepath;
        if strcmp(paramsField, 'exp_filepath')
            loadExpFile(filepath);
        elseif strcmp(paramsField, 'cal_filepath')
            loadCalFile(filepath);
        end
    end

    function selectFolder(paramsField)
        % Set main figure invisible to prevent file dialog being hidden
        fig.Visible = 'off';
        disp('    SELECT FOLDER TO SAVE RESULTS IN...');
        folderpath = uigetdir(pwd, 'Select folder to save results in');
        fig.Visible = 'on';
        if isequal(folderpath, 0), return; end  % User canceled selection
        folderpath = fullfile(folderpath);
        params.(paramsField) = folderpath;
        TextFieldSave.Value = folderpath;
        TextFieldSave.Tooltip = folderpath;
        % Mark setup flag 3 as completed
        params.setupAllCompleted(3) = 1;
        fprintf('Analysis results will be saved to: %s\n', folderpath);
    end

    function on_Close(fig)
        selection = uiconfirm(fig,"Closing this window will abort the script! Continue?",...
            "Confirmation");
        switch selection
            case 'OK'
                delete(fig);
            case 'Cancel'
                return
        end
    end


end