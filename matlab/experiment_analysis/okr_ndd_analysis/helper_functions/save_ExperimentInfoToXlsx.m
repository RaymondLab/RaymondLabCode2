function save_ExperimentInfoToXlsx(params, nblocks)
%SAVE_EXPERIMENTINFOTOXLSX Summary of this function goes here
%   Detailed explanation goes here

% Get calibration scaling factors
cal_customscaling = 'No';
if ~isempty(params.customScaleChs)
    cal_customscaling = 'Yes';
    cal_scalechs = sprintf('(%g, %g)', params.customScaleChs(1), params.customScaleChs(2));
elseif ~isempty(params.calScaleChs)
    cal_scalechs = sprintf('(%g, %g)', params.calScaleChs(1), params.calScaleChs(2));
else
    cal_scalechs = '(none, none)';
end

% Get the username 
username = {getenv('USERNAME'), getenv('USER'), 'None'};
username = username{find(~cellfun('isempty', username), 1)};

% Get the network hostname
[~, net_hostname] = system('hostname');
net_hostname = strtrim(net_hostname);   % remove trailing newline

% Get experiment date
if isfield(params,'exp_date') && ~isempty(params.exp_date)
    exp_date = params.exp_date;
else
    exp_date = '';
end

KEY = {'analysis'; 'analysis_date'; ...
    'exp_filepath'; 'exp_date'; 'num_blocks'; 'train_type'; ...
    'exp_subid'; 'exp_subcond'; 'exp_subcohort'; 'exp_taskcond'; ...
    'stim_frequency'; 'sampling_rate'; 'lp_cutoff'; 'sg_window'; ...
    'cal_filepath'; 'cal_scalechs'; 'cal_customscaling'; ...
    'desaccade_method'; 'saccade_threshold'; 'saccade_lrpad'; 'min_goodlen'; ...
    'username'; 'comp_serialnum'; 'comp_os'; 'net_hostname'; 'matlab_ver'};

VALUE = {'ndd_okrAnalysis'; string(datetime('now')); ...
    params.exp_filepath; exp_date; nblocks; params.exp_traintype; ...
    params.exp_subid; params.exp_subcond; params.exp_subcohort; params.exp_taskcond;
    params.exp_stimfreq; params.fs; params.lowpassCutoff; params.filterWindow; ...
    params.cal_filepath; cal_scalechs; cal_customscaling; ...
    params.saccadeMethod; params.saccadeThresh; params.saccadeLRPad; params.minGoodChunk_len; ...
    username; getSerialNumber(); computer; net_hostname; version};

xlsx = table(KEY, VALUE);

% Write this data to new sheets in excel sheet
[folder,file,~] = fileparts(params.save_filepath);
xlsx_filepath = fullfile(folder, [file, '.xlsx']);
writetable(xlsx, xlsx_filepath, 'Sheet','ExpInfo');

    %% Helper functions
    function serial = getSerialNumber()
        if ispc
            [~, result] = system('wmic bios get serialnumber');
            serial = strtrim(regexp(result, '(?<=SerialNumber\s*)\S+', 'match', 'once'));
        elseif ismac
            [~, result] = system('ioreg -rd1 -c IOPlatformExpertDevice | awk -F''"'' ''/IOPlatformSerialNumber/{print $4}''');
            serial = strtrim(result);
        elseif isunix
            [~, serial] = system('cat /sys/class/dmi/id/product_serial 2>/dev/null');
            serial = strtrim(serial);
        else
            serial = '';
        end
    
        if isempty(serial)
            serial = 'Unknown';
        end
    end

end