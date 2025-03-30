function [vars] = manualDesaccade_APP(vars, mode, position)
%MANUALDESACCADE_APP

%% Load data
loadAnalysisInfo_APP;

if ~isstring(position)
    %% Set line interval as a saccade
    % Get the indices of the corresponding left and right times of the line
    t = vid.time_upsampled_aligned;
    [~, leftIdx] = min(abs(t - position(1,1)));
    [~, rightIdx] = min(abs(t - position(2,1)));
    
    % Set the selected interval as a saccade
    switch mode
        case "mag1"
            mag1.saccades_manual(leftIdx:rightIdx) = 1;
        case "mag2"
            mag2.saccades_manual(leftIdx:rightIdx) = 1;
        case "vid"
            vid.saccades_manual(leftIdx:rightIdx) = 1;
        otherwise
            error('"'+mode+'" is an invalid mode argument for manualDesaccade_APP');
    end
else
    % Reset the manually added saccades
    switch mode
        case "mag1"
            mag1.saccades_manual(1:end) = 0;
        case "mag2"
            mag2.saccades_manual(1:end) = 0;
        case "vid"
            vid.saccades_manual(1:end) = 0;
        otherwise
            error('"'+mode+'" is an invalid mode argument for manualDesaccade_APP');
    end
end

%% Save data
saveAnalysisInfo_APP;

end