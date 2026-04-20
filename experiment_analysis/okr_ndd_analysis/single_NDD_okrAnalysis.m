clear; clc; close all;
params.script_folderpath = fileparts(which(matlab.desktop.editor.getActiveFilename));
addpath(genpath(params.script_folderpath));

% OKR_NDD_EXP_ANALYSIS - OKR Normal-Dim-Dimmer (NDD) Experiment Analysis
% Important: This script requires "RaymondLabCode2" and its subfolders added to MATLAB path!


%% SCRIPT PARAMETERS (Do not edit unless you know what you are doing!)

% General experiment parameters
params.exp_stimfreq  = 1.0;    % Training stimulus frequency (Default 1.0)
params.exp_traintype = 'OKR';  % Training type (Default 'OKR')
% 
% Define timepoint blocks and corresponding group ids
% Do not adjust this for now as the group_NDD_okrAnalysis is not set up for this
params.timepoint_ids    = [2, 3, 4, 14, 15, 16, 59, 60, 61];  % Default [2, 3, 4, 59, 60, 61]
params.timepoint_groups = [0, 0, 0, 15, 15, 15, 60, 60, 60];  % Default [0, 0, 0, 60, 60, 60]

% Subject condition options for dropdown menu
params.subcond_options = {'None', 'NR', '23DMR', '24DMR', 'DR'}; 

% Task (training) condition options for dropdown menu
params.taskcond_options = {'None', ...
    'TNOR', 'TDIM', 'TDMR', ...
    '2TNOR', '2TDIM', '2TDMR', ...
    'D0', 'D20', 'D40', 'D60', ...
    'D60N20', 'D60N40', ...
    'D120', 'D120N20'}; 

% Number of perturbations to run for F-Tests
params.ftestPerturbs = 10000;  % (Default 10000, or 5000 for tests)

% Custom channel scaling is used if tuple is provided (Default [])
params.customScaleChs = [];

% Data preprocessing parameters
params.filtering_options = {'Butter', 'FIF'};
params.transientThresh   = 1;   % Transient "spike" removal threshold (Default 1)
params.lowpassCutoff     = 15;  % Lowpass filter cutoff frequency (Default 11)
params.filterWindow      = 11;  % SG filter window size (Default 30)

% FIF-specific settings
params.fifSettings.delta = 0.1;  
params.fifSettings.Xi = 2;  
params.fifSettings.NIMFs = 4;  
params.fifSettings.alpha = 24;
params.fifSettings.minIMF = 4;

% Saccade removal parameters
params.saccadeMethod    = 'MAD';  % Saccade detection method (Default 'MAD' or 'SVT')
params.saccadeThresh    = 7;   % Saccade detection threshold (Default 7 or 1000)
params.saccadeLRPad     = 0.05;   % Padding in seconds to expand saccade candidates (Default 0.05)
params.minGoodChunk_len = 200;     % Minimum allowed gap between saccades (Default 200)

% File and folder paths (by default, all should be '')
params.exp_filepath    = '';  % Experiment file (.smr or .mat) path
params.cal_filepath    = '';  % Calibration file (.mat) path
params.save_folderpath = '';  % Target save folder path

% LPC x saccade-threshold comparison sweep (opt-in: needs >=2 values in BOTH)
params.lowpass_cutoffs = [7, 15, 30];   % e.g. [9, 15, 30]; empty or scalar = no comparison
params.saccade_threshs = [2.5, 5, 9];   % e.g. [2.5, 5, 9];  empty or scalar = no comparison
params.ntColor = [1 0 0];      % NT peak marker color (robustness tabs)
params.tnColor = [0 0 1];      % TN peak marker color (robustness tabs)


%% BEGINNING OF SCRIPT SETUP

disp('--- RUNNING OKR NDD ANALYSIS ! ---');

% Run initial analysis setup GUI
[params, data] = run_okrNDD_setupGUI(params);

% If GUI was closed prematurely and not completed, throw an error to abort script
if ~all(params.setupAllCompleted)
    error('--- Analysis window was manually closed. ABORTING SCRIPT! ---');
else
    fprintf('\nAnalysis setup was successfully completed!\n')
end

% Extract the data from the GUI
nBlocks = data.nBlocks;
blockIds = data.blockIds;
blockTypes = data.blockTypes;
bssTimes = data.bssTimes;
cycleLength = data.cycleLength;
cycleTimes = data.cycleTimes;
chairvel_raw = data.chairvel_raw;
drumvel_raw = data.drumvel_raw;
hevel_raw = data.hevel_raw;
hepos_raw = data.hepos_raw;
saccadePad = data.saccadePad;
% NOTE: `data` is intentionally kept in scope so the opt-in lpc/gc
% comparison sweep (see bottom of script) can reuse it. It is cleared
% later, after all analyses complete.


%% BLOCK-LEVEL ANALYSIS

tic;
fprintf('\nStarting analysis:\n');
disp('    Analyzing data by blocks...');
tpbn = 0;  % Counter for keeping track of timepoint blocks
clear blocks;
for ii = 1:nBlocks
    bk = run_blockAnalysis(params, data, ii, params.lowpassCutoff, params.saccadeThresh);
    if ismember(ii, params.timepoint_ids)
        tpbn = tpbn + 1;
        bk.timePoint = sprintf('T%d', params.timepoint_groups(tpbn));
    else
        bk.timePoint = 'NA';
    end
    blocks(ii) = bk;
end

% Clean up workspace
clearvars('-except', 'analysis', 'blocks', 'cycleLength', 'cycleTimes', 'data', 'diffdata', 'params', 'timepoints');

% Generate Excel sheet of general block analysis results
if ~isequal(params.save_folderpath, 0) & isfield(params, 'save_filepath')
    disp('    Saving block data results to Excel sheet...');
    save_BlockAnalysisToXlsx(blocks, params);
end


%% TIMEPOINT GROUP ANALYSIS
disp('    Analyzing data grouped by timepoint...');
timepoints = run_timepointAnalysis(blocks, params, cycleLength);
ntp = length(timepoints);

% Clean up workspace
clearvars('-except', 'analysis', 'blocks', 'cycleLength', 'cycleTimes', 'data', 'diffdata', 'ntp', 'params', 'timepoints');


%% DIFFDATA ANALYSIS
disp('    Conducting diffdata analysis...');
diffdata = run_diffdataAnalysis(timepoints);

% Add sheets to existing Excel for diffdata analysis results
if ~isequal(params.save_folderpath, 0) & isfield(params, 'save_filepath')
    disp('    Adding diffdata analysis results to Excel file...');
    save_DiffDataAnalysisToXlsx(timepoints, diffdata, params);
    disp('    Adding experiment information metadata to Excel file...');
    save_ExperimentInfoToXlsx(params, length(blocks));
end

% Clean up workspace
clearvars('-except', 'analysis', 'blocks', 'data', 'diffdata', 'params', 'timepoints');

disp('    Compiling analysis results...');
analysis.params     = params;
analysis.blocks     = blocks;
analysis.timepoints = timepoints;
analysis.diffdata   = diffdata;
clear blocks timepoints diffdata;

fprintf('\nGenerating plots of analysis results:\n');
fignum = plot_standard_subplots(analysis, params.save_folderpath);
fignum = plot_okr_NDD_analysis(analysis, params.save_folderpath, fignum);

% Opt-in lowpass-cutoff x saccade-threshold comparison sweep
if numel(params.lowpass_cutoffs) >= 2 && numel(params.saccade_threshs) >= 2
    fprintf('    Running LPC x saccade-threshold comparison sweep...\n');
    [T_lpcgc, D_lpcgc] = run_lpcgcComparison(params, data);
    fignum = plot_lpcgcComparison(analysis, T_lpcgc, D_lpcgc, params.save_folderpath, fignum);
    clear T_lpcgc D_lpcgc;
end

temp_analysis = analysis;
if ~isequal(params.save_folderpath, 0) & isfield(params, 'save_filepath')
    fprintf('    Saving analysis parameters to: %s\n', params.save_filepath); 
    clear analysis;
    analysis.params = temp_analysis.params;
    save(params.save_filepath, 'analysis');
    clear analysis;
    analysis = temp_analysis;
end

% Final cleanup
clearvars('-except', 'analysis');

fprintf('\n---ANALYSIS COMPLETED!---\n\n');
msgbox('Analysis complete!', 'Done');
toc
