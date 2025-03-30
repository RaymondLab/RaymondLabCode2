function gainAnalysisApp()
%GAINANALYSISAPP Summary of this function goes here
%   Detailed explanation goes here

% If window already exists, change focus to window and exit early
persistent fig
if ~isempty(fig) && ishandle(fig) && isvalid(fig)
    figure(fig);
    return
end

% Get resolution of smallest available display
resolution_scale = 0.75;
monitors = get(groot,'MonitorPositions');
[~, minIdx] = min(prod(monitors(:,3:end), 2));
maxwh = monitors(minIdx, 3:4);
figwh = round(maxwh * resolution_scale);
figxy = round((maxwh - figwh) * 0.5);
uixywh = [figxy, figwh];
clear resolution_scale monitors minIdx figwh figxy

%% Create instance of data model
F = Spike2_Model();

%% Initialize new launcher figure
fig = uifigure('Toolbar','none', ... 
    'MenuBar','none', ...
    'NumberTitle','off', ...
    'Name','Gain Analysis App', ...
    'Position',uixywh, ...
    'WindowState','maximized');

% Preallocate ui structs to hold figure components
ui1 = struct;
ui2 = struct;
ui3 = struct;
analysis = struct;
btable = [];

%% Set up Window 1
ui1.GridLayout = uigridlayout(fig, ...
    'RowHeight',{105,50,'1x'}, ...
    'ColumnWidth',{'1x'}, ...
    'Padding',[20 20 20 20]);
% Provide access to main figure "fig" by adding reference in UserData
ui1.GridLayout.UserData = fig;
% Separate custom components are created per grid in GridLayout
% Add component for loading smr recording and calibration files
ui1.LoadFiles = LoadSmrCalFiles_Comp(ui1.GridLayout);
% Add compoenent for adjusting Butterworth and Savitzky-Golay filter parameters
ui1.FilterSettings = FilterSettings_Comp(ui1.GridLayout);
ui1.FilterSettings.setLowpassCutoff(F.lowpassCutoff);
ui1.FilterSettings.setSavgolWindow(F.savgolWindow);
% Add component to display plots of filtered data
ui1.FilteredDataPlots = FilteredDataPlots_Comp(ui1.GridLayout);
% Group listeners into cell array for easy cleanup
% Listeners will automatically call functions when triggered by the given event
ui1.Listeners = { ...
    addlistener(ui1.LoadFiles, "smrLoaded", @(src,evt) F.setSmr(src.smr)), ...
    addlistener(ui1.LoadFiles, "scaleChLoaded", @(src,evt) F.setScaleCh(src.scaleCh)), ...
    addlistener(ui1.LoadFiles, "scaleChLoaded", @(src,evt) ui1.FilterSettings.updateDropdown(F.nblocks)), ...
    addlistener(ui1.LoadFiles, "nextPushed", @(src,evt) loadBlocksWindow()), ...
    addlistener(ui1.FilterSettings, "lowpassCutoffChanged", @(src,evt) F.setLowpassCutoff(src.lowpassCutoff)), ...
    addlistener(ui1.FilterSettings, "savgolWindowChanged", @(src,evt) F.setSavgolWindow(src.savgolWindow)), ...
    addlistener(ui1.FilterSettings, "blockNumChanged", @(src,evt) ui1.FilteredDataPlots.updateAll(F.data, F.blocktimes, F.samplerate, ui1.FilterSettings.blockNum)), ...
    addlistener(F, "dataChanged", @(src,evt) ui1.FilteredDataPlots.updateAll(F.data, F.blocktimes, F.samplerate, ui1.FilterSettings.blockNum)), ...
    addlistener(F, "filterValueChanged", @(src,evt) ui1.FilteredDataPlots.updateFiltered(F.data.hepos, F.data.hevel, F.blocktimes, F.samplerate, ui1.FilterSettings.blockNum))
};


%% Subfunctions  
    % Set up Window 2
    function loadBlocksWindow()
        % Display progress dialog
        pdlg = uiprogressdlg(fig, ...
            'Title','Loading, please wait...');
        % Reset figure and clear ui1 to free memory
        ui1.GridLayout.delete;
        clear ui1;
        % Initialize new window
        ui2.GridLayout = uigridlayout(fig, ...
            'RowHeight',{105,'1x'}, ...
            'ColumnWidth',{'1x'}, ...
            'Padding',[20 20 20 20]);
        % Provide access to main figure "fig" by adding reference in UserData
        ui2.GridLayout.UserData = fig;
        % Add components to display plots of block data
        ui2.BlockSettingsPanel = BlockSettingsPanel_Comp(ui2.GridLayout);
        ui2.BlockSettingsPanel.updateBlockSettings(F.schemaidx, F.nblocks, F.blocktypes);
        ui2.BlockDataPlots = BlockDataPlots_Comp(ui2.GridLayout);
        ui2.BlockDataPlots.updateTraces(F.data);
        ui2.BlockDataPlots.updateBlockRegions(F.blocktimes, F.blocktypes);
        ui2.Listeners = addlistener(ui2.BlockSettingsPanel, ...
            "nextPushed", ...
            @(src,evt) loadDesaccadeWindow());
        close(pdlg);
    end

    function loadDesaccadeWindow()
        % Display progress dialog
        pdlg = uiprogressdlg(fig, ...
            'Title','Loading, please wait...');
        % Reset figure and clear ui2 to free memory
        btable = ui2.BlockDataPlots.blockTable;
        F.setStimulusFrequency(btable.Frequency(1));
        ui2.GridLayout.delete;
        clear ui2;
        % Initialize new window
        ui3.GridLayout = uigridlayout(fig, ...
            'RowHeight',{105,105,'1x'}, ...
            'ColumnWidth',{'1x'}, ...
            'Padding',[20 20 20 20]);
        % Provide access to main figure "fig" by adding reference in UserData
        ui3.GridLayout.UserData = fig;
        % Add components to display plots of saccade data
        ui3.SaccadesSettingsPanel = SaccadesSettingsPanel_Comp(ui3.GridLayout);
        ui3.AnalysisSettingsPanel = AnalysisSettingsPanel_Comp(ui3.GridLayout);
        ui3.Listeners = addlistener(ui3.SaccadesSettingsPanel, "saccadesPushed", @(src,evt) plotReviewSaccades());
        close(pdlg);
    end

    function plotReviewSaccades()
        % Saccade detection method parameters
        methodIdx = ui3.SaccadesSettingsPanel.DropDown.ValueIndex;
        methodThresh = ui3.SaccadesSettingsPanel.NumericFields(1).Value;
        saccWin = ui3.SaccadesSettingsPanel.NumericFields(2).Value;
        minLen = ui3.SaccadesSettingsPanel.NumericFields(3).Value;
        samplerate = F.samplerate;
        stimFreq = F.stimulusFrequency;
        D = F.data;
        
        % Apply corresponding saccade method
        if isequal(methodIdx, 1)     % Moving Median MAD
            [saccades, ~, hevel, ~] = desaccadeMAD(D.heposraw, samplerate, stimFreq, saccWin, saccWin, methodThresh, minLen);
        elseif isequal(methodIdx, 2) % Squared Velocity Threshold
            [saccades, ~, hevel, ~] = desaccadeSVT(D.heposraw, samplerate, stimFreq, saccWin, saccWin, methodThresh, minLen);
        end

        % % Loop through each block
        % for ii = 1:F.nblocks
        %     bstart = round(F.blocktimes(ii,1) * F.samplerate);
        %     bend = round(F.blocktimes(ii,2) * F.samplerate);
        % 
        %     % Extract block data
        %     chairVel = D.hhvel(bstart:bend);
        %     drumVel = D.htvel(bstart:bend);
        %     eyeVel = hevel(bstart:bend);
        % 
        %     blockLength = length(eyeVel);
        %     blockTime = (1:blockLength) / samplerate;
        % 
        %     % Calculate fit coefficients
        %     y1 = sin(2 * pi * stimFreq * blockTime(:));
        %     y2 = cos(2 * pi * stimFreq * blockTime(:));
        %     constant = ones(blockLength, 1);
        %     vars = [y1 y2 constant];
        % 
        %     b = regress(chairVel, vars);
        %     chairVel_amp = sqrt(b(1)^2 + b(2)^2);
        %     chairVel_angle = rad2deg(atan2(b(2),b(1)));
        % 
        %     b = regress(drumVel, vars);
        %     drumVel_amp = sqrt(b(1)^2 + b(2)^2);
        %     drumVel_angle = rad2deg(atan2(b(2),b(1)));
        % 
        %     b = regress(eyeVel(~saccades), vars(~saccades,:));
        %     eyeVel_amp = sqrt(b(1)^2 + b(2)^2);
        %     eyeVel_phase = rad2deg(atan2(b(2),b(1)));
        % 
        %     % Calculate eye gain as VOR or OKR based on head signal
        %     % Chair/VOR Stimulus
        %     if chairVel_amp > 3
        %         reference_amp = chairVel_amp;
        %         reference_angle = chairVel_angle;
        %         idealEyeVel = drumVel - chairVel;
        %     % Drum/OKR Stimulus
        %     elseif drumVel_amp > 3
        %         reference_amp = drumVel_amp;
        %         reference_angle = drumVel_angle;
        %         idealEyeVel = drumVel;
        %     % No Motor Stimulus
        %     else
        %         reference_amp = 1;
        %         reference_angle = 0;
        %         idealEyeVel = zeros(1, blockLength);
        %     end
        % 
        %     % Eye calculations relative to drum/chair
        %     eyeVel_rel_gain = eyeVel_amp / reference_amp;
        %     eyeVel_rel_phase = (eyeVel_phase - reference_angle);
        %     eyeVel_rel_phase = mod(eyeVel_rel_phase,360) - 180;
        %     eyeVel_cycleFit = sin(2*pi*freq*cycleTime + deg2rad(eyeVel_rel_phase+180)) * eyeVel_amp;
        % 
        %     % Calculate Averages and Metrics
        %     startpt = max(1, round(mod(-reference_angle,360)/360 * samplerate/stimFreq));
        % 
        %     [eyeVel_des_mat, eyeVel_des_cycleMean] = VOR_breakTrace(cycleLength, startpt, eyeVel_des);
        %     [~, headVel_cycleMean]                 = VOR_breakTrace(cycleLength, startpt, headVel);
        %     [~, drumVel_cycleMean]                 = VOR_breakTrace(cycleLength, startpt, drumVel);
        %     [omit_mat, ~]                          = VOR_breakTrace(cycleLength, startpt, double(omitCenters));
        %     [~, idealEye_cycleMean]                = VOR_breakTrace(cycleLength, startpt, idealEyeVel);
        % 
        %     % Calculate Extras
        %     badCycles       = any(omit_mat,2);
        %     goodCount       = sum(~badCycles);
        %     if size(badCycles,1)==0
        %         goodCount = 1 ;
        %     end
        %     eyeVel_des_Sem  = nanstd(eyeVel_des_mat)./sqrt(sum(~isnan(eyeVel_des_mat)));
        % 
        %     if goodCount > 0
        %         eyeVel_good_cycleMean = nanmean(eyeVel_des_mat(~badCycles,:), 1);
        %         eyeVel_good_cycleStd  = nanstd(eyeVel_des_mat(~badCycles,:));
        %     else
        %         eyeVel_good_cycleMean = zeros(size(eyeVel_des_mat(1,:)));
        %         eyeVel_good_cycleStd  = zeros(size(eyeVel_des_mat(1,:)));
        %     end
        % end
        % 
        % % Save results and information in the analysis struct
        % analysis.saccades = saccades;
    end
    
    % function closeRequestFunc(src, event)
    % % 'CloseRequestFcn',@closeRequestFunc, ...
    % selection = uiconfirm(src,'Save current changes?', 'Confirm Close'); 
    % switch selection 
    %     case 'OK'
    %         delete(src)
    %     case 'Cancel'
    %         return 
    % end
    % end

end