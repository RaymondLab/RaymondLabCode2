function gainAnalysisApp()
%GAINANALYSISAPP Summary of this function goes here
%   Detailed explanation goes here

% If window already exists, change focus to window and exit early
persistent fig
if ~isempty(fig) && ishandle(fig) && isvalid(fig)
    figure(fig);
    return
end

% Get resolution of smallest available display
resolution_scale = 0.75;
monitors = get(groot,'MonitorPositions');
[~, minIdx] = min(prod(monitors(:,3:end), 2));
maxwh = monitors(minIdx, 3:4);
figwh = round(maxwh * resolution_scale);
figxy = round((maxwh - figwh) * 0.5);
uixywh = [figxy, figwh];
clear resolution_scale monitors minIdx figwh figxy

%% Initialize new launcher figure
fig = uifigure('Toolbar','none', ... 
    'MenuBar','none', ...
    'NumberTitle','off', ...
    'Name','Gain Analysis App', ...
    'Position',uixywh);

%% Create tab group container for app tabs
tabgp = uitabgroup(fig, 'Units','normalized', 'Position',[0 0 1 1]);

%% Create tabs for each step of analysis
tab1.uitab = uitab(tabgp, 'Title','Load');
tab2.uitab = uitab(tabgp, 'Title','Process');
tab3.uitab = uitab(tabgp, 'Title','Clean');
tab4.uitab = uitab(tabgp, 'Title','Analyze');

%% Initialize layout of "Load" tab (tab 1)
tab1.grid.uigrid = uigridlayout(tab1.uitab, ...
    'RowHeight',{105,'1x'}, ...
    'ColumnWidth',{'1x'});

% Grid square for "Load" files components
tab1.grid.gridFolder.uigrid = uigridlayout(tab1.grid.uigrid, ...
    'RowHeight',{'1x','1x'}, ...
    'ColumnWidth',{100,'1x'});
tab1.grid.gridFolder.uibutton{1} = uibutton(tab1.grid.gridFolder.uigrid, ...
    'Text','Load Recording', ...
    'ButtonPushedFcn',@(src,evt) changeFilepath(1,1));
tab1.grid.gridFolder.uieditfield{1} = uieditfield(tab1.grid.gridFolder.uigrid, ...
    'InputType','text', ...
    'Placeholder',' Path to Spike2 recording file to analyze', ...
    'ValueChangedFcn',@(src,evt) changeFilepath(2,1));
tab1.grid.gridFolder.uibutton{2} = uibutton(tab1.grid.gridFolder.uigrid, ...
    'Text','Load Calibration', ...
    'ButtonPushedFcn',@(src,evt) changeFilepath(1,2));
tab1.grid.gridFolder.uieditfield{2} = uieditfield(tab1.grid.gridFolder.uigrid, ...
    'InputType','text', ...
    'Placeholder',' Path to calibration .mat file containing scaling factors to use', ...
    'ValueChangedFcn',@(src,evt) changeFilepath(2,2));

% Grid square for displaying drum, chair, and eye position channels
tab1.grid.gridRawPlots.uigrid = uigridlayout(tab1.grid.uigrid, ...
    'RowHeight',{'1x','1x','1x'}, ...
    'ColumnWidth',{'1x'});

tab1.grid.gridRawPlots.uiax{1} = uiaxes(tab1.grid.gridRawPlots.uigrid);
tab1.grid.gridRawPlots.uiax{2} = uiaxes(tab1.grid.gridRawPlots.uigrid);
tab1.grid.gridRawPlots.uiax{3} = uiaxes(tab1.grid.gridRawPlots.uigrid);
tab1.grid.gridRawPlots.uiax{4} = uiaxes(tab1.grid.gridRawPlots.uigrid);


%% App variables
% Initialize struct to store app component-related data
% This struct gets reset each time the app is run
app = struct;

% Initialize struct to store parameter-related data
% Load previous "vars" data if it already exists in the target folder
vars = struct;

% Initialize struct to store recording data
data = struct;


%% Subfunctions
    function changeFilepath(fieldtype, num)
        if isequal(fieldtype, 1)
            extensions = {{'.smr'}, {'.mat'}};
            fig.Visible = 'off';
            [filename,folderpath] = uigetfile(extensions{num}, 'Select Spike2 recording file', cd);
            fig.Visible = 'on';
            % Early exit if filepath is zero (i.e. uigetfile was cancelled)
            if isequal(filename, 0)
                return
            end
            filepath = fullfile(folderpath, filename);
        elseif isequal(fieldtype, 2)
            filepath = tab1.grid.gridFolder.uieditfield{num}.Value;
        else
            error('Invalid fieldtype passed to changeFilePath function');
        end

        % Throw error if filepath is not a valid file
        if ~isfile(filepath)
            error('Invalid file path: %s', filepath);
        end

        if isequal(num, 1) && ~isempty(tab1.grid.gridFolder.uieditfield{num}.Value)
            % Provide confirmation message of resetting state of app
            selection = uiconfirm(fig, ...
                'Loading this recording file will reset the app. Continue?', ...
                'Confirm Load New File', ...
                'Icon','warning');
            % Do nothing and exit early if cancelled
            if selection == "Cancel"
                return
            end
            % Reset the app ui whenever recording file is changed
            resetAppUI();
        elseif isequal(num, 2)
            % Verify calibration file contains valid scaling factors
            cal = load(filepath, 'scaleCh1', 'scaleCh2');
            % Verify calibration file has both scaling factors
            status1 = isfield(cal,'scaleCh1') && isfield(cal,'scaleCh2');
            try
                % Verify only one scaling factor is nonzero
                status2 = isequal(nnz([cal.scaleCh1,cal.scaleCh2]), 1);
            catch
                status2 = 0;
            end
            if ~status1 | ~status2
                error('Calibration .mat file has invalid scaling factors: %s', filepath);
            end
            vars.scaleCh1 = cal.scaleCh1;
            vars.scaleCh2 = cal.scaleCh2;
        end
        
        % Set file and folder information
        [app.folders{num},vars.filename{num},vars.extensions{num}] = fileparts(filepath);
        tab1.grid.gridFolder.uieditfield{num}.Value = filepath;

        % Update current folder and import data if loading recording file
        if isequal(num, 1)
            cd(app.folders{num});
            if isfile([vars.filename{num},'.mat'])
                data = load([vars.filename{num},'.mat'], 'smr');
            else
                data.smr = importSpike2File(filepath);
            end
            loadDataInfo();
        end

        % Final verification that everything is good before proceeding
        good = [~isempty(tab1.grid.gridFolder.uieditfield{1}.Value), ...
            ~isempty(tab1.grid.gridFolder.uieditfield{2}.Value)];
        if all(good)
            %todo
        end
    end

    function loadDataInfo()
        ch_htvel = find(ismember(data.smr.channeltitles, 'htvel'));        
        if isequal(length(ch_htvel), 1)
            plotChannelDataUI(tab1.grid.gridRawPlots.uiax{1}, data.smr, ch_htvel);
        end

        ch_hhvel = find(ismember(data.smr.channeltitles, 'hhvel'));
        if isequal(length(ch_hhvel), 1)
            plotChannelDataUI(tab1.grid.gridRawPlots.uiax{2}, data.smr, ch_hhvel);
        end

        ch_hepos1 = find(strcmpi(data.smr.channeltitles, 'hepos1'));
        if isequal(length(ch_hepos1), 1)
            plotChannelDataUI(tab1.grid.gridRawPlots.uiax{3}, data.smr, ch_hepos1);
        end

        ch_hepos2 = find(strcmpi(data.smr.channeltitles, 'hepos2'));
        if isequal(length(ch_hepos2), 1)
            plotChannelDataUI(tab1.grid.gridRawPlots.uiax{4}, data.smr, ch_hepos2);
        end
    end

    function resetAppUI()
        % Reset app, vars, and data structs
        app = struct;
        app.folders = cell(2, 1);
        vars = struct;
        vars.filenames = cell(2, 1);
        vars.extensions = cell(2, 1);
        vars.scaleCh1 = nan;
        vars.scaleCh2 = nan;
        data = struct;

        % Reset app components
        tab1.grid.gridFolder.uieditfield{1}.Value = '';
        tab1.grid.gridFolder.uieditfield{2}.Value = '';
    end

    % function closeRequestFunc(src, event)
    % % 'CloseRequestFcn',@closeRequestFunc, ...
    % selection = uiconfirm(src,'Save current changes?', 'Confirm Close'); 
    % switch selection 
    %     case 'OK'
    %         delete(src)
    %     case 'Cancel'
    %         return 
    % end
    % end

end