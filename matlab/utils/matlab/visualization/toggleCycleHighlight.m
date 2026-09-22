function toggleCycleHighlight(evt, hLines)
    clicked = evt.Peer;
    if ~any(clicked == hLines)
        % Non-cycle-line entry: preserve default toggle-visibility behavior
        clicked.Visible = matlab.lang.OnOffSwitchState(~strcmp(clicked.Visible,'on'));
        return;
    end
    wasHighlighted = strcmp(clicked.Tag, 'highlighted');
    set(hLines, 'LineWidth', 1, 'Tag', '');   % reset all three
    if ~wasHighlighted
        set(clicked, 'LineWidth', 3.5, 'Tag', 'highlighted');
    end
end