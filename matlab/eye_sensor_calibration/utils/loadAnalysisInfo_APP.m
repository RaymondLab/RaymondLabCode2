%% load analysis information
try
    % don't load from a file if you don't have to
    mag1 = vars.mag1;
    mag2 = vars.mag2;
    vid = vars.vid;
    if ~isfield(mag1, 'align_sign'); mag1.align_sign = 1.0; end
    if ~isfield(mag2, 'align_sign'); mag2.align_sign = 1.0; end
    if ~isfield(vid,  'align_sign'); vid.align_sign  = 1.0; end
    if ~isfield(vid, 'stim') || ~isfield(vid, 'stim_label') || isempty(vid.stim)
        [~, filenameroot]= fileparts(cd);
        fullfilename = fullfile(cd, [filenameroot,'.smr']);
        stimdata = importSpike(fullfilename, {'hhvel', 'htvel', 'Keyboard'});
        lightpulses = stimdata(end).data;
        segmentStart = lightpulses(1);
        segmentEnd = segmentStart+vid.time(end);
        stimdata = resettime(datseg(stimdata, [segmentStart segmentEnd]));
        hhvel = double(stimdata(1).data);
        [hhvel_amp,~,~,~,~] = fit_sineWave(hhvel, mag1.samplerate, vars.stimFreq);
        htvel = double(stimdata(2).data);
        [htvel_amp,~,~,~,~] = fit_sineWave(htvel, mag1.samplerate, vars.stimFreq);
        if hhvel_amp > 3
            vid.stim = hhvel;
            vid.stim_label = 'hhvel Channel';
        elseif htvel_amp > 3
            vid.stim = htvel;
            vid.stim_label = 'htvel Channel';
        else
            vid.stim = nan(size(hhvel));
            vid.stim_label = '';
        end
    end
catch
    try
        pathname = cd;
        [~, filenameroot]= fileparts(pathname);
        load(fullfile(cd, [filenameroot '_analysis.mat']))
    catch
        error('Unable to read desaccading data')
    end
end