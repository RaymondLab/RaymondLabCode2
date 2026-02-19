function save_BlockAnalysisToXlsx(blocks, params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
cols = {'blockType', 'blockNumber','stimFreq', 'startTime', 'endTime', 'hevel_rel_gain', 'hevel_rel_phase', 'headamp', 'headangle', ...
    'drumamp', 'drumangle', 'hevel_amp', 'hevel_phase', 'saccadeFrac', 'rsquare', 'pval', 'variance', 'nGoodCycles', ...
    'good_rel_gain', 'good_rel_gain_SEM', 'good_rel_phase', 'good_rel_phase_SEM', 'good_amp', 'good_amp_SEM', 'good_phase', 'good_phase_SEM', ...
    'stimvel_amp', 'stimvel_amp_SEM', 'stimvel_phase', 'stimvel_phase_SEM'};
T = struct2table(blocks);
T = T(:,cols);
T = renamevars(T, ["blockType", "blockNumber", "stimFreq", "startTime", "endTime", "hevel_rel_gain", "hevel_rel_phase", "hevel_amp", "hevel_phase"], ...
    ["Type", "Timepoint", "Frequency", "StartTime", "EndTime", "eyeHgain", "eyeHphase", "eyeHamp", "eyeHangle"]);
T.Timepoint = T.Timepoint - 1;
T = addvars(T, NaN(size(T.Timepoint)), 'NewVariableNames','Notes', 'After','EndTime');

% To add a gap between the "Notes" and "eyeHgain" columns, we split the table
gapAfterCol = 6;
T1 = T(:, 1:gapAfterCol);
T2 = T(:, gapAfterCol+1:end);

% Write each table to the Excel file with a gap between them
[folder,file,~] = fileparts(params.save_filepath);
xlsx_savepath = fullfile(folder, [file,'.xlsx']);  
startCol = char('A' + gapAfterCol + 3);  % +1 for the gap
if isfile(xlsx_savepath), delete(xlsx_savepath); end
writetable(T1, xlsx_savepath, 'Sheet','Blocks', 'Range','A1');
writetable(T2, xlsx_savepath, 'Sheet','Blocks', 'Range',[startCol,'1']);

end