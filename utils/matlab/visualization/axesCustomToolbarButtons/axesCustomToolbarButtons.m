function axesCustomToolbarButtons(ax, interactions)
%AXESCUSTOMTOOLBARBUTTONS Replaces axes toolbar with custom set of buttons
%   Replaces the default axes toolbar with a custom set of more useful
%   buttons for exploring plotted time series data.
%
%   Description of the provided buttons:
%   "Restore View"    - This is a default button that restores the original 
%       view of axes or tiled chart layouts.
%   "Horizontal Zoom" - Toggles horizontal-only zooming.
%   "Vertical Zoom"   - Toggles vertical-only zooming.
%
% Inputs:
%   ax  : Target axes to apply this custom toolbar
%
%   ----------------------------------------------------------------------
%   Author: Brian Angeles, Stanford University, 01/2025
%   ----------------------------------------------------------------------

% Manually set default interactions of axes
if ~exist('interactions', 'var')
    interactions = [panInteraction, zoomInteraction];
end
ax.Interactions = interactions;

% Replace axes toolbar with only the "restor view" button
tb = axtoolbar(ax, {'restoreview'});

% Add a "horizontal zoom" state button
hBtn = axtoolbarbtn(tb, 'state', 'Icon','tool_zoom_x.gif');
hBtn.ValueChangedFcn = {@horizontalZoomState, ax};
hBtn.Tooltip = 'Toggle Horizontal-Only Zoom Mode';

% Add a "vertical zoom" state button
vBtn = axtoolbarbtn(tb, 'state', 'Icon','tool_zoom_y.gif');
vBtn.ValueChangedFcn = {@verticalZoomState, ax};
VBtn.Tooltip = 'Toggle Vertical-Only Zoom Mode';

% Subfunctions that define functionality of custom buttons
    function horizontalZoomState(src, ~, ax)
        switch src.Value
            case 'off'
                ax.Interactions = [panInteraction, ...
                    zoomInteraction('Dimensions','xy')];
            case 'on'
                vBtn.Value = 'off';
                ax.Interactions = [panInteraction, ...
                    zoomInteraction('Dimensions','x')];
        end
    end

    function verticalZoomState(src, ~, ax)
        switch src.Value
            case 'off'
                ax.Interactions = [panInteraction, ...
                    zoomInteraction('Dimensions','xy')];
            case 'on'
                hBtn.Value = 'off';
                ax.Interactions = [panInteraction, ...
                    zoomInteraction('Dimensions','y')];
        end
    end
end