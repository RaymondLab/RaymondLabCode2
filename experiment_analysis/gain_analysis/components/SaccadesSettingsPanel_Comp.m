classdef SaccadesSettingsPanel_Comp < matlab.ui.componentcontainer.ComponentContainer
    %SACCADESSETTINGSPANEL_COMP Summary of this class goes here
    %   Detailed explanation goes here

    % Define public events that can be used to trigger outside functions
    events
        saccadesPushed
    end

    properties
        DropDown matlab.ui.control.DropDown
        NumericFields (1,:) matlab.ui.control.NumericEditField
    end

    properties (Access = private)
        Panel matlab.ui.container.Panel
        GridLayout matlab.ui.container.GridLayout
        DropDownLabel matlab.ui.control.Label
        FieldLabels (1,:) matlab.ui.control.Label
        SaccadesButton matlab.ui.control.Button
    end

    methods (Access = protected)
        function setup(obj)
            % Define Panel and GridLayout inside container
            obj.Panel = uipanel(obj.Parent);
            obj.Panel.Layout.Row = 1;
            gridColumnWidth = {'fit',225,'0.5x','fit',40,'0.5x','fit',40,'0.5x','fit',45,'1x',200};
            obj.GridLayout = uigridlayout(obj.Panel, ...
            'RowHeight',{'1x',25,'1x'}, ...
            'ColumnWidth',gridColumnWidth, ...
            'Padding',[30 5 5 5]);

            obj.DropDownLabel = uilabel(obj.GridLayout, ...
                'Text','Desaccade Method: ', ...
                'FontSize',15, ...
                'HorizontalAlignment','right');
            obj.DropDownLabel.Layout.Row = 2;
            obj.DropDownLabel.Layout.Column = 1;
            obj.DropDown = uidropdown(obj.GridLayout, ...
                'FontSize',15, ...
                'Items',["Moving Median MAD", "Squared Velocity Threshold"]);
            obj.DropDown.Layout.Row = 2;
            obj.DropDown.Layout.Column = 2;
            
            obj.FieldLabels(1) = uilabel(obj.GridLayout, ...
                'Text','Threshold: ', ...
                'FontSize',15, ...
                'HorizontalAlignment','right');
            obj.FieldLabels(1).Layout.Row = 2;
            obj.FieldLabels(1).Layout.Column = 4;
            obj.NumericFields(1) = uieditfield(obj.GridLayout, ...
                "numeric", ...
                'FontSize',15, ...
                'Value',7);
            obj.NumericFields(1).Layout.Row = 2;
            obj.NumericFields(1).Layout.Column = 5;

            obj.FieldLabels(2) = uilabel(obj.GridLayout, ...
                'Text','Saccade Window (ms): ', ...
                'FontSize',15, ...
                'HorizontalAlignment','right');
            obj.FieldLabels(2).Layout.Row = 2;
            obj.FieldLabels(2).Layout.Column = 7;
            obj.NumericFields(2) = uieditfield(obj.GridLayout, ...
                "numeric", ...
                'FontSize',15, ...
                'Value',50);
            obj.NumericFields(2).Layout.Row = 2;
            obj.NumericFields(2).Layout.Column = 8;

            obj.FieldLabels(3) = uilabel(obj.GridLayout, ...
                'Text','Min. Good Length (ms): ', ...
                'FontSize',15, ...
                'HorizontalAlignment','right');
            obj.FieldLabels(3).Layout.Row = 2;
            obj.FieldLabels(3).Layout.Column = 10;
            obj.NumericFields(3) = uieditfield(obj.GridLayout, ...
                "numeric", ...
                'FontSize',15, ...
                'Value',500);
            obj.NumericFields(3).Layout.Row = 2;
            obj.NumericFields(3).Layout.Column = 11;
           

            % Create Review Saccades button
            obj.SaccadesButton = uibutton(obj.GridLayout, ...
                'Text',{'Review','Saccades'}, ...
                'FontSize',18, ...
                'Enable','on', ...
                'ButtonPushedFcn',@(src,evt) obj.reviewSaccadesPushed());
            obj.SaccadesButton.Layout.Column = length(gridColumnWidth);
            obj.SaccadesButton.Layout.Row = [1 3];
        end

        function update(~)
            % Required by ComponentContainer
            % Do nothing
        end
    end

    methods (Access = private)
        function reviewSaccadesPushed(obj)
            % Notify listeners that the "next" button was Pushed
            notify(obj, 'saccadesPushed');
        end 
    end

end