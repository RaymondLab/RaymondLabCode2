classdef AnalysisSettingsPanel_Comp < matlab.ui.componentcontainer.ComponentContainer
    %ANALYSISSETTINGSPANEL_COMP Summary of this class goes here
    %   Detailed explanation goes here

    % Define public events that can be used to trigger outside functions
    events
        startPushed
    end

    properties (Access = private)
        Panel matlab.ui.container.Panel
        GridLayout matlab.ui.container.GridLayout
        DropDownLabel matlab.ui.control.Label
        DropDown matlab.ui.control.DropDown
        FieldLabels (1,:) matlab.ui.control.Label
        NumericFields (1,:) matlab.ui.control.NumericEditField
        StartButton matlab.ui.control.Button
    end

    methods (Access = protected)
        function setup(obj)
            % Define Panel and GridLayout inside container
            obj.Panel = uipanel(obj.Parent);
            obj.Panel.Layout.Row = 2;
            gridColumnWidth = {'1x','fit',300,'1x',200};
            obj.GridLayout = uigridlayout(obj.Panel, ...
            'RowHeight',{'1x',25,'1x'}, ...
            'ColumnWidth',gridColumnWidth, ...
            'Padding',[20 5 5 5]);

            % Create dropdown component
            obj.DropDownLabel = uilabel(obj.GridLayout, ...
                'Text','Analysis/Project: ', ...
                'FontSize',15, ...
                'HorizontalAlignment','right');
            obj.DropDownLabel.Layout.Row = 2;
            obj.DropDownLabel.Layout.Column = 2;
            obj.DropDown = uidropdown(obj.GridLayout, ...
                'FontSize',15, ...
                'Items',["Default (Sine)", "Default (Step)"]);
            obj.DropDown.Layout.Row = 2;
            obj.DropDown.Layout.Column = 3;

            % Create Review Saccades button
            obj.StartButton = uibutton(obj.GridLayout, ...
                'Text','START', ...
                'FontSize',18, ...
                'Enable','on', ...
                'ButtonPushedFcn',@(src,evt) obj.startButtonPushed());
            obj.StartButton.Layout.Column = length(gridColumnWidth);
            obj.StartButton.Layout.Row = [1 3];

            obj.Panel.Enable = 'off';
        end

        function update(~)
            % Required by ComponentContainer
            % Do nothing
        end
    end

    methods
        % function updateBlockSettings(obj, schemaidx, nblocks, btypes)
        %     obj.FieldLabels(1).Text = sprintf('Block Schema:  %d', schemaidx);
        %     obj.FieldLabels(2).Text = sprintf('Total number of blocks:  %d', nblocks);
        %     obj.FieldLabels(3).Text = sprintf('Number of pre/post-test blocks:  %d', sum(btypes==1));
        %     obj.FieldLabels(4).Text = sprintf('Number of training blocks:  %d', sum(btypes==2));
        %     obj.FieldLabels(5).Text = sprintf('Number of mid-test blocks:  %d', sum(btypes==4));
        % end
    end

    methods (Access = private)
        function startButtonPushed(obj)
            % Notify listeners that the "next" button was Pushed
            notify(obj, 'saccadesPushed');
        end 
    end

end