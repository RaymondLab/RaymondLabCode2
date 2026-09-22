classdef BlockSettingsPanel_Comp < matlab.ui.componentcontainer.ComponentContainer
    %BLOCKSETTINGSPANEL_COMP Summary of this class goes here
    %   Detailed explanation goes here

    % Define public events that can be used to trigger outside functions
    events
        nextPushed
    end

    properties (Access = private)
        Panel matlab.ui.container.Panel
        GridLayout matlab.ui.container.GridLayout
        FieldLabels (1,:) matlab.ui.control.Label
        NextButton matlab.ui.control.Button
    end

    methods (Access = protected)
        function setup(obj)
            % Define Panel and GridLayout inside container
            obj.Panel = uipanel(obj.Parent);
            obj.Panel.Layout.Row = 1;
            gridColumnWidth = {'1x',200,'1x',200,'1x',200,'1x',200,'1x',200,'1x',200};
            obj.GridLayout = uigridlayout(obj.Panel, ...
            'RowHeight',{'1x',40,'1x'}, ...
            'ColumnWidth',gridColumnWidth);

            obj.FieldLabels(1) = uilabel(obj.GridLayout, ...
                'Text','Block Schema:  0', ...
                'FontSize',17, ...
                'FontWeight','bold', ...
                'HorizontalAlignment','center');
            obj.FieldLabels(1).Layout.Row = 2;
            obj.FieldLabels(1).Layout.Column = 2;

            obj.FieldLabels(2) = uilabel(obj.GridLayout, ...
                'Text','Total blocks:  0', ...
                'FontSize',16, ...
                'HorizontalAlignment','center');
            obj.FieldLabels(2).Layout.Row = 2;
            obj.FieldLabels(2).Layout.Column = 4;

            obj.FieldLabels(3) = uilabel(obj.GridLayout, ...
                'Text','Pre/Post-test blocks:  0', ...
                'FontSize',16, ...
                'FontColor',[0 0.6 0], ...
                'HorizontalAlignment','center');
            obj.FieldLabels(3).Layout.Row = 2;
            obj.FieldLabels(3).Layout.Column = 6;

            obj.FieldLabels(4) = uilabel(obj.GridLayout, ...
                'Text','Training blocks:  0', ...
                'FontSize',16, ...
                'FontColor',[0.6 0 0], ...
                'HorizontalAlignment','center');
            obj.FieldLabels(4).Layout.Row = 2;
            obj.FieldLabels(4).Layout.Column = 8;

            obj.FieldLabels(5) = uilabel(obj.GridLayout, ...
                'Text','Mid-test blocks:  0', ...
                'FontSize',16, ...
                'FontColor',[0 0 0.6], ...
                'HorizontalAlignment','center');
            obj.FieldLabels(5).Layout.Row = 2;
            obj.FieldLabels(5).Layout.Column = 10;

            % Create next button
            obj.NextButton = uibutton(obj.GridLayout, ...
                'Text','NEXT', ...
                'FontWeight','bold', ...
                'FontSize',18, ...
                'Enable','on', ...
                'ButtonPushedFcn',@(src,evt) obj.nextButtonPushed());
            obj.NextButton.Layout.Row = [1 3];
            obj.NextButton.Layout.Column = length(gridColumnWidth);
        end

        function update(~)
            % Required by ComponentContainer
            % Do nothing
        end
    end

    methods
        function updateBlockSettings(obj, schemaidx, nblocks, btypes)
            obj.FieldLabels(1).Text = sprintf('Block Schema:  %d', schemaidx);
            obj.FieldLabels(2).Text = sprintf('Total blocks:  %d', nblocks);
            obj.FieldLabels(3).Text = sprintf('Pre/Post-test blocks:  %d', sum(btypes==1));
            obj.FieldLabels(4).Text = sprintf('Training blocks:  %d', sum(btypes==2));
            obj.FieldLabels(5).Text = sprintf('Mid-test blocks:  %d', sum(btypes==4));
        end
    end

    methods (Access = private)
        function nextButtonPushed(obj)
            % Notify listeners that the "next" button was Pushed
            notify(obj, 'nextPushed');
        end        
    end

end