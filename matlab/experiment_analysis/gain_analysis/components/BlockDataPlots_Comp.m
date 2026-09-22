classdef BlockDataPlots_Comp < matlab.ui.componentcontainer.ComponentContainer
    %BLOCKDATAPLOTS_COMP Summary of this class goes here
    %   Detailed explanation goes here

    properties (SetAccess = private)
        blockTable (:,:) table = table('Size',[1 3], ...
            'VariableNames',{'Block', 'Label', 'Frequency'}, ...
            'VariableTypes',{'double', 'string', 'double'});
        blabels (:,1) cell = {"VORD"; "X2"; "VORD"; "VORD"};
    end

    properties (Access = private, Transient, NonCopyable)
        GridLayout matlab.ui.container.GridLayout
        axs (1,:) matlab.ui.control.UIAxes
        axlines (1,:) matlab.graphics.chart.primitive.Line
        axregions (1,:) cell
        axtable matlab.ui.control.Table
    end

    methods (Access = protected)
        function setup(obj) 
            % Define GridLayout inside container
            obj.GridLayout = uigridlayout(obj.Parent, ...
                'RowHeight',{'1x','1x','1x'}, ...
                'ColumnWidth',{'1x',180});
            obj.GridLayout.Layout.Row = 2;

            obj.axtable = uitable(obj.GridLayout, ...
                'Data',obj.blockTable,...
                'ColumnEditable',[false true true], ...
                'ColumnWidth','fit', ...
                'DisplayDataChangedFcn',@obj.tableDataChanged);
            obj.axtable.Layout.Row = [1 3];
            obj.axtable.Layout.Column = 2;
            addStyle(obj.axtable, uistyle('HorizontalAlignment','left'));

            % Initialize uiaxes for velocity traces
            for ii = 1:3
                obj.axs(ii) = uiaxes(obj.GridLayout);
                obj.axs(ii).Layout.Row = ii;
                obj.axs(ii).Layout.Column = 1;
                obj.axs(ii).Interactions = [panInteraction, zoomInteraction];
                axesCustomToolbarButtons(obj.axs(ii));
                obj.axlines(ii) = plot(obj.axs(ii), NaN, NaN, ...
                    'Color',[0 0 0 0.4], 'LineWidth',1, 'HandleVisibility','off');
                switch ii
                    case 1
                        ylabel(obj.axs(ii), 'downsampled htvel (deg/s)', 'FontSize',14);
                    case 2
                        ylabel(obj.axs(ii), 'downsampled hhvel (deg/s)', 'FontSize',14);
                    case 3
                        xlabel(obj.axs(ii), 'Time (s)', 'FontSize',14); 
                        ylabel(obj.axs(ii), 'downsampled hevel (deg/s)', 'FontSize',14);
                end
                box(obj.axs(ii), 'on');
                grid(obj.axs(ii), 'on');
                obj.axs(ii).GridAlpha = 0.07;
            end
            linkaxes(obj.axs, 'x');
        end

        function update(~)
            % Required by ComponentContainer
            % Do nothing
        end
    end

    methods 
        function updateTraces(obj, data)
            % Downsample data to improve performance
            dsf = 5;
            % Time x-axis data
            times = downsample(seconds(data.Time), dsf);
            % Update velocity traces
            obj.axlines(1).XData = times;
            obj.axlines(1).YData = downsample(data.HTVEL, dsf);
            obj.axlines(2).XData = times;
            obj.axlines(2).YData = downsample(data.HHVEL, dsf);
            obj.axlines(3).XData = times;
            obj.axlines(3).YData = downsample(data.hevel, dsf);
            % Update axes limits
            xlim(obj.axs(3), [times(1), times(end)]);
            ylim(obj.axs(1), [-25 25]);
            ylim(obj.axs(2), [-25 25]);
            ylim(obj.axs(3), [-60 60]);
        end

        function updateBlockRegions(obj, btimes, btypes)
            bcolors = 'grgb';
            bcolors2 = {[0 0.6 0], [1 0 0], [0 0.6 0], [0 0 1]};
            nblocks = length(btypes);
            btable = table('Size',[nblocks 3], ...
                'VariableNames',{'Block','Label','Frequency'}, ...
                'VariableTypes',{'double','string','double'});
            btable.Block(:) = 1:nblocks;
            btable.Label(:) = obj.blabels(btypes);
            btable.Frequency(:) = ones(nblocks, 1);
            obj.UserData = btypes;
            obj.blockTable = btable;
            obj.axtable.Data = btable;
            for ii = 1:4
                maskii = (btypes == ii);
                if gt(sum(maskii), 0)
                    [r,c] = find(maskii);
                    addStyle(obj.axtable, ...
                        uistyle('FontColor',bcolors2{ii}), ...
                        'cell',[r,c+1]);
                    obj.axregions{ii} = [];
                    obj.axregions{ii} = xregion(obj.axs(1), ...
                        btimes(maskii,1), ...
                        btimes(maskii,2), ...
                        'FaceColor',bcolors(ii), ...
                        'FaceAlpha',0.4);
                end
            end
        end
    end

    methods (Access = private)
        function tableDataChanged(obj, ~, evt)
            btypes = obj.UserData; 
            evtRow = evt.DisplaySelection(1);
            evtCol = evt.DisplaySelection(2);
            newValue = obj.axtable.Data{evtRow,evtCol};
            switch evtCol
                case 2
                    typeIdx = btypes(evtRow);
                    obj.blabels{typeIdx} = newValue;
                    obj.blockTable.Label(:) = obj.blabels(btypes);
                case 3
                    obj.blockTable.Frequency(:) = zeros(length(btypes),1) + newValue;
            end
            obj.axtable.Data = obj.blockTable;
        end
    end
    
end