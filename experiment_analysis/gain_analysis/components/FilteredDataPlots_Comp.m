classdef FilteredDataPlots_Comp < matlab.ui.componentcontainer.ComponentContainer
    %FILTEREDDATAPLOTS_COMP Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = private, Transient, NonCopyable)
        GridLayout matlab.ui.container.GridLayout
        axs (1,:) matlab.ui.control.UIAxes
        axlines (1,:) matlab.graphics.chart.primitive.Line
    end

    methods (Access = protected)
        function setup(obj) 
            % Use brewermap to get three diverging colors
            c = brewermap(10, 'Seaborn');
            c1 = [c(8,:), 0.2]; 
            c2 = [c(1,:), 0.9];
            c3 = [c(2,:), 0.8];
            % Define GridLayout inside container
            obj.GridLayout = uigridlayout(obj.Parent, ...
                'RowHeight',{'1x'}, ...
                'ColumnWidth',{'1x','1x'});
            obj.GridLayout.Layout.Row = 3;
            % Initialize uiaxes for eye position traces
            obj.axs(1) = uiaxes(obj.GridLayout);
            axesCustomToolbarButtons(obj.axs(1));
            obj.axlines(1) = plot(obj.axs(1), NaN, NaN, ...
                'Color',c1, 'LineWidth',3, ...
                'DisplayName',' raw eye position');
            hold(obj.axs(1), 'on');
            obj.axlines(2) = plot(obj.axs(1), NaN, NaN, ...
                'Color',c2, 'LineWidth',2, ...
                'DisplayName',' 40Hz-filtered eye position');
            obj.axlines(3) = plot(obj.axs(1), NaN, NaN, ...
                'Color',c3, 'LineWidth',2.5, ...
                'DisplayName',' filtered eye position');
            hold(obj.axs(1), 'off');
            xlabel(obj.axs(1), 'Time (s)', 'FontSize',14);
            ylabel(obj.axs(1), 'hepos (deg)', 'FontSize',14);
            box(obj.axs(1), 'on');
            grid(obj.axs(1), 'on');
            obj.axs(1).GridAlpha = 0.07;
            legend(obj.axs(1), 'Location','northwest', 'FontSize',12);
            % Initialize uiaxes for eye velocity traces
            obj.axs(2) = uiaxes(obj.GridLayout);
            axesCustomToolbarButtons(obj.axs(2));
            obj.axlines(4) = plot(obj.axs(2), NaN, NaN, ...
                'Color',c1, 'LineWidth',3, ...
                'DisplayName',' raw eye velocity');
            hold(obj.axs(2), 'on');
            obj.axlines(5) = plot(obj.axs(2), NaN, NaN, ...
                'Color',c2, 'LineWidth',2.5, ...
                'DisplayName',' 40Hz-filtered, 3ms-windowed eye velocity'); 
            obj.axlines(6) = plot(obj.axs(2), NaN, NaN, ...
                'Color',c3, 'LineWidth',2.5, ...
                'DisplayName',' filtered eye velocity');
            hold(obj.axs(2), 'off');
            xlabel(obj.axs(2), 'Time (s)', 'FontSize',14);
            ylabel(obj.axs(2), 'hevel (deg/s)', 'FontSize',14);
            box(obj.axs(2), 'on');
            grid(obj.axs(2), 'on');
            obj.axs(2).GridAlpha = 0.07;
            legend(obj.axs(2), 'Location','northwest', 'FontSize',12);
            linkaxes([obj.axs(1), obj.axs(2)], 'x');
        end

        function update(~)
            % Required by ComponentContainer
            % Do nothing
        end
    end

    methods 
        function updateAll(obj, data, btimes, fs, bnum)
            % Time x-axis data
            times = seconds(data.Time);
            % Generate array mask for given block number
            % Invalid block numbers will simply plot all the data
            try
                if lt(bnum, 1)
                    maskids = 1:length(times);
                else
                    maskids = round(btimes(bnum,:)*fs);
                    maskids = maskids(1,1):maskids(1,2);
                    times = times(maskids);
                end
            catch
                maskids = 1:length(times);
            end
            % Update position traces
            obj.axlines(1).XData = times;
            obj.axlines(1).YData = data.heposraw(maskids);
            obj.axlines(2).XData = times;
            obj.axlines(2).YData = data.heposdefault(maskids);
            obj.axlines(3).XData = times;
            obj.axlines(3).YData = data.hepos(maskids);
            % Update velocity traces
            obj.axlines(4).XData = times;
            obj.axlines(4).YData = data.hevelraw(maskids);
            obj.axlines(5).XData = times;
            obj.axlines(5).YData = data.heveldefault(maskids);
            obj.axlines(6).XData = times;
            obj.axlines(6).YData = data.hevel(maskids);
            xlim(obj.axs(2), [times(1), times(end)]);
            ylim(obj.axs(2), [-300, 300]);
        end

        function updateFiltered(obj, hepos, hevel, btimes, fs, bnum)
            if lt(bnum, 1)
                % Invalid block numbers will attempt to plot all the data
                obj.axlines(3).YData = hepos;
                obj.axlines(6).YData = hevel;
            else
                % Generate array mask for given block number
                maskids = round(btimes(bnum,:)*fs);
                maskids = maskids(1,1):maskids(1,2);
                obj.axlines(3).YData = hepos(maskids);
                obj.axlines(6).YData = hevel(maskids);
            end
        end
    end
end