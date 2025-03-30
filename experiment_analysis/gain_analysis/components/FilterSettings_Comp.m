classdef FilterSettings_Comp < matlab.ui.componentcontainer.ComponentContainer
    %FILTERSETTINGS_COMP Summary of this class goes here
    %   Detailed explanation goes here

    % Define public events that can be used to trigger outside functions
    events
        blockNumChanged
        lowpassCutoffChanged
        savgolWindowChanged
    end

    properties (SetAccess = protected)
        blockNum (1,:) char {mustBeText} = 1;
        % Cutoff frequency (Hz) used for lowpass Butterworth filter
        lowpassCutoff (1,1) double {mustBeNumeric,mustBePositive} = 11.0;
        % Window (ms) used for movingslopeCausal ("savgol") filter
        savgolWindow (1,1) double {mustBeInteger,mustBePositive} = 30;
        % Flag to bypass listeners when updating observed properties
        isUpdating (1,1) logical = false;
    end

    properties (Access = private)
        Panel matlab.ui.container.Panel
        GridLayout matlab.ui.container.GridLayout
        BlockNumLabel matlab.ui.control.Label
        BlockNumDropDown matlab.ui.control.DropDown
        LowpassLabel matlab.ui.control.Label
        LowpassNumericField matlab.ui.control.NumericEditField
        SavgolLabel matlab.ui.control.Label
        SavgolNumericField matlab.ui.control.NumericEditField
    end

    methods (Access = protected)
        function setup(obj)
            % Define Panel and GridLayout inside container
            obj.Panel = uipanel(obj.Parent);
            obj.Panel.Layout.Row = 2;
            obj.GridLayout = uigridlayout(obj.Panel, ...
            'RowHeight',{'1x'}, ...
            'ColumnWidth',{'1x',225,55,'0.5x',225,55,'1x',225,55,'1x'});
            % Lowpass cutoff label and numeric field components
            obj.LowpassLabel = uilabel(obj.GridLayout, ...
                'Text','Lowpass Cutoff Frequency (Hz): ', ...
                'FontSize',14, ...
                'HorizontalAlignment','right');
            obj.LowpassLabel.Layout.Column = 2;
            obj.LowpassNumericField = uieditfield(obj.GridLayout, ...
                "numeric", ...
                'FontSize',14, ...
                'Value',obj.lowpassCutoff, ...
                'ValueChangedFcn',@obj.cutoffChanged);
            % Block number label and dropdown components
            obj.BlockNumLabel = uilabel(obj.GridLayout, ...
                'Text','Block Number: ', ...
                'FontSize',14, ...
                'HorizontalAlignment','right');
            obj.BlockNumLabel.Layout.Column = 5;
            obj.BlockNumDropDown = uidropdown(obj.GridLayout, ...
                'Items',{'all', '1'}, ...
                'ValueIndex',2, ...
                'FontSize',14, ...
                'ValueChangedFcn',@obj.dropdownChanged);
            % Sav-Gol window label and numeric field components
            obj.SavgolLabel = uilabel(obj.GridLayout, ...
                'Text','Differentiation Window Size (ms): ', ...
                'FontSize',14, ...
                'HorizontalAlignment','right');
            obj.SavgolLabel.Layout.Column = 8;
            obj.SavgolNumericField = uieditfield(obj.GridLayout, ...
                "numeric", ...
                'FontSize',14, ...
                'Value',obj.savgolWindow, ...
                'ValueChangedFcn',@obj.windowChanged);
        end

        function update(~)
            % Required by ComponentContainer
            % Do nothing
        end
    end

    methods
        function updateDropdown(obj, nblocks)
            if gt(nblocks, 0)
                newItems = ['all', arrayfun(@(x) num2str(x),1:nblocks,'UniformOutput',false)];
            else 
                newItems = {'all'};
                obj.blockNum = 0;
            end
            obj.BlockNumDropDown.Items = newItems;
        end

        function setLowpassCutoff(obj, value)
            if obj.isUpdating || lt(value, 0.1) || gt(value, 499.9)
                return; % Exit early if new value is invalid
            end
            obj.isUpdating = true;
            obj.LowpassNumericField.Value = value;
            obj.lowpassCutoff = value;
            obj.isUpdating = false;
        end

        function setSavgolWindow(obj, value)
            if obj.isUpdating || lt(value, 3) || gt(value, 9999)
                return; % Exit early if new value is invalid
            end
            obj.isUpdating = true;
            obj.SavgolNumericField.Value = value;
            obj.savgolWindow = value;
            obj.isUpdating = false;
        end
    end
    
    methods (Access = private)
        function cutoffChanged(obj, ~, evt)
            newValue = evt.Value;
            if obj.isUpdating || lt(newValue, 0.1) || gt(newValue, 499.9)
                return; % Exit early if new value is invalid
            end
            % Set new lowpass cutoff value
            obj.lowpassCutoff = newValue;
            % Notify listeners that the lowpass cutoff was changed
            notify(obj, 'lowpassCutoffChanged');
        end

        function dropdownChanged(obj, ~, evt)
            newValue = round(evt.ValueIndex - 1);
            if obj.isUpdating
                return; % Exit early if flag is on
            end
            % Set new block number value
            obj.blockNum = newValue;
            % Notify listeners that the lowpass cutoff was changed
            notify(obj, 'blockNumChanged');
        end

        function windowChanged(obj, ~, evt)
            newValue = round(evt.Value); % Ensure new value is an integer
            if obj.isUpdating || lt(newValue, 3) || gt(newValue, 9999)
                return; % Exit early if new value is invalid
            end
            % Set new savgol window value
            obj.savgolWindow = newValue;
            % Notify listeners that the savgol window was changed
            notify(obj, 'savgolWindowChanged');
        end
    end

end