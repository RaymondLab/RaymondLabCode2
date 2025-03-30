classdef LoadSmrCalFiles_Comp < matlab.ui.componentcontainer.ComponentContainer
    %LOADSMRCALFILES_COMP Summary of this class goes here
    %   Detailed explanation goes here

    % Define public events that can be used to trigger outside functions
    events
        smrLoaded
        scaleChLoaded
        nextPushed
    end

    properties (SetAccess = protected)
        % Spike2 struct file generated via `importSpike2File` function
        smr (1,1) struct
        % Eye calibration scaling factors of respective magnet channels
        scaleCh (2,1) double {mustBeNumeric} = [NaN; NaN];
    end

    properties (Access = private)
        GridLayout matlab.ui.container.GridLayout
        LoadSmrButton matlab.ui.control.Button
        LoadCalButton matlab.ui.control.Button
        SmrTextField matlab.ui.control.EditField
        CalTextField matlab.ui.control.EditField
        NextButton matlab.ui.control.Button
    end

    methods (Access = protected)
        function setup(obj)
            % Define GridLayout inside container
            obj.GridLayout = uigridlayout(obj.Parent, ...
                'RowHeight',{'1x','1x'}, ...
                'ColumnWidth',{130,'1x',200});
            obj.GridLayout.Layout.Row = 1;
            % Create button and text field for loading recording file
            obj.LoadSmrButton = uibutton(obj.GridLayout, ...
                'Text','Load Recording', ...
                'FontSize',14, ...
                'ButtonPushedFcn',@(src,evt) obj.selectFile('.smr'));
            obj.SmrTextField = uieditfield(obj.GridLayout, ...
                'InputType','text', ...
                'Placeholder',' Load Spike2 recording .smr file to analyze', ...
                'FontSize',14, ...
                'Editable','off');
            % Create button and text field for loading calibration file
            obj.LoadCalButton = uibutton(obj.GridLayout, ...
                'Text','Load Calibration', ...
                'FontSize',14, ...
                'ButtonPushedFcn',@(src,evt) obj.selectFile('*cal*.mat'));
            obj.LoadCalButton.Layout.Row = 2;
            obj.LoadCalButton.Layout.Column = 1;
            obj.CalTextField = uieditfield(obj.GridLayout, ...
                'InputType','text', ...
                'Placeholder',' Load calibration .mat file containing scaling factors', ...
                'FontSize',14, ...
                'Editable','off');
            obj.CalTextField.Layout.Row = 2;
            obj.CalTextField.Layout.Column = 2;
            % Create next button
            obj.NextButton = uibutton(obj.GridLayout, ...
                'Text','NEXT', ...
                'FontWeight','bold', ...
                'FontSize',18, ...
                'Enable','off', ...
                'ButtonPushedFcn',@(src,evt) obj.nextButtonPushed());
            obj.NextButton.Layout.Row = [1 2];
            obj.NextButton.Layout.Column = 3;
        end

        function update(~)
            % Required by ComponentContainer
            % Do nothing
        end
    end

    methods (Access = private)
        function selectFile(obj, fileExt)
            % Set main figure invisible to prevent file dialog being hidden
            obj.Parent.UserData.Visible = 'off';
            [filename,folderpath] = uigetfile(fileExt, ['Select a ',fileExt,' file'], cd);
            obj.Parent.UserData.Visible = 'on';
            if isequal(filename, 0)
                return; % User canceled selection
            end
            filepath = fullfile(folderpath, filename);
            % Display progress dialog
            pdlg = uiprogressdlg(obj.Parent.UserData, ...
                'Title','Loading, please wait...');
            % Load respective file
            if fileExt == ".smr"
                obj.loadSmrFile(filepath)
            elseif fileExt == "*cal*.mat"
                obj.loadCalFile(filepath)
            end
            % Close progress dialog
            close(pdlg);
        end
        
        function loadSmrFile(obj, filepath)
            try
                % Give confirmation message that changing smr file resets the app
                if ~isempty(fieldnames(obj.smr))
                    selection = uiconfirm(obj.Parent.UserData, ...
                        'Loading a new .smr file will reset the app. Continue?', ...
                        'Confirm Load New File', ...
                        'Icon','warning');
                    if selection == "Cancel"
                        return; % User canceled selection
                    end
                end
                [folderpath,filename,~] = fileparts(filepath);
                cd(folderpath);
                if isfile([filename,'.mat'])
                    % Check if a .mat version of the file already exists
                    data = load([filename,'.mat']).smr;
                else
                    % Otherwise import from the .smr version
                    data = importSpike2File(filepath);
                end
                % Verify imported data is valid
                if ~isfield(data, 'filename')
                    badfilename = [filename,'.mat'];
                    warning('Folder contains invalid smr .mat file: %s', badfilename);
                    return;
                end
                % Set text field and smr property
                obj.SmrTextField.Value = filepath;
                obj.smr = data;
                % Enable "next" button if both text fields are not empty
                obj.enableNextButton();
                % Notify listeners that the smr property was loaded
                notify(obj, 'smrLoaded');
            catch ME
                warning(ME.identifier, 'Error loading file: %s', ME.message);
            end
        end

        function loadCalFile(obj, filepath)
            try
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
                    warning('Calibration .mat file has invalid scaling factors: %s', filepath);
                    return;
                end
                % Set text field and scaleCh property
                obj.CalTextField.Value = filepath;
                obj.scaleCh = [cal.scaleCh1; cal.scaleCh2];
                % Enable "next" button if both text fields are not empty
                obj.enableNextButton();
                % Notify listeners that the scaleCh property was loaded
                notify(obj, 'scaleChLoaded');
            catch ME
                warning(ME.identifier, 'Error loading file: %s', ME.message);
            end
        end

        function enableNextButton(obj)
            if ~isempty(obj.SmrTextField.Value) & ~isempty(obj.CalTextField.Value)
                obj.NextButton.Enable = 'on';
            else
                obj.NextButton.Enable = 'off';
            end
        end

        function nextButtonPushed(obj)
            % Notify listeners that the "next" button was Pushed
            notify(obj, 'nextPushed');
        end
    end

end