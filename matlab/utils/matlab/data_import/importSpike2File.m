function smr = importSpike2File(fullfilepath, savefile)
%IMPORTSPIKE2FILE Imports Spike2 file as a MATLAB struct
%   Imports a Spike2 file into an smr struct variable that is structured 
%   for maximum compatibility and performance in both MATLAB and Python. 
%   This function is still a work inprogress and should be updated to
%   resolve any bugs or errors may occur.
% 
%   Requires the MATLAB SON32 library version 2.4 or higher.
% 
% Inputs:
%   filepath : string representing full path to Spike2 file.
%   savefile : boolean that determines whether to save smr as a .mat file.
%
% Outputs:
%   smr      : structure that contains imported data from Spike2 file.
%      .filename       : original filename of Spike2 file.
%      .fileextension  : original file extension (i.e. smr or smrx).
%      .date           : date and time of recording.
%      .filecomments   : 5-by-1 cell of 5 comments saved in file.
%      .channelnumbers : N-by-1 array of N valid channel numbers.
%      .channelkinds   : N-by-1 array of N valid channel kinds.
%      .channeltitles  : N-by-1 cell of N valid channel titles.
%      .channels       : N-by-1 cell of N valid channel data.
%
%   Vector channel data are S-by-1 arrays of S values.
%   Non-vector channel data are given as S-by-M arrays of M sets of S values.
%
%   smr is saved as a .mat file that can be loaded into Python using Scipy:
%   --- Example Python code snippet
%   import scipy.io as sio
%   smr = sio.loadmat('example.mat', squeeze_me=True, struct_as_record=False)
%   ---
%   See the Scipy documentation for more information:
%   https://docs.scipy.org/doc/scipy/tutorial/io.html#matlab-files
%
%   ----------------------------------------------------------------------
%   Author: Brian Angeles, Stanford University, 01/2025
%   ----------------------------------------------------------------------

% Display open file dialog box if filepath not given
if ~exist('fullfilepath','var') || isempty(fullfilepath) || ~isfile(fullfilepath)
    [filename, filepath] = uigetfile('*.smr');
    fullfilepath         = fullfile(filepath, filename);
end

% Save a .mat file by default
if ~exist('savefile', 'var')
    savefile = true;
end

% Initialize output smr variable as a struct
smr = struct;

% Extract filename and extension from the filepath
[filepath, smr.filename, smr.fileextension] = fileparts(fullfilepath);

% Ensure file is a valid Spike2 recording file
if strcmpi(smr.fileextension, '.smr') == 1
    fid = fopen(fullfilepath, 'r', 'l');
else
    error('Error. %s is not a valid Spike2 file\n', filename);
end

if fid < 0
    error('Error. File not found, please check path');
else
    fprintf('Importing file: %s%s\n', smr.filename, smr.fileextension);
end

% Get file header and channel information
H = SONFileHeader(fid);
C = SONChanList(fid);

% Preallocate and set initial file properties
smr.date = '';
smr.filecomments   = H.fileComment;
smr.channelnumbers = [];
smr.channelkinds   = [];
smr.channeltitles  = {};

% Get channel number, kind, and port arrays
channelnumbers = [C.number]';
channelkinds   = [C.kind]';
channelports   = [C.phyChan]';

% Read the Spike2 file and get the number of valid channels
F         = readSpikeFile(fullfilepath, channelnumbers);
nchannels = length(F);
channels  = cell(nchannels, 1);
for ii = 1:nchannels
    % Get data and header of ii-th channel
    fdata            = F(ii).data;
    fheader          = F(ii).header;
    % Preallocate and set channel properties
    chan.number      = fheader.channel;
    chanIdx          = find(channelnumbers==chan.number, 1);
    chan.kind        = channelkinds(chanIdx);
    chan.port        = channelports(chanIdx);
    chan.title       = fheader.title;
    chan.comment     = fheader.comment;
    chan.units       = '';
    chan.samplerate  = NaN;
    chan.data        = [];
    chan.markertimes = [];
    chan.markerdata  = {};
    % Not sure if we still have files with these channel types
    % If so, we can modify the fillGaps function to support them
    if strcmp(fheader.channeltype, 'Episodic Waveform')
        error('Unsupported channel type: Episodic Waveform');
        %d(iIn) = fillGaps(d(iIn));
    end
    switch chan.kind
        case {1, 9}  % Waveform (ADC) or RealWave channels
            chan.units = fheader.adc.Units;
            chan.samplerate = round(1/prod(fheader.adc.SampleInterval), 1);
            chan.data = (single(fdata.adc) * fheader.adc.Scale) + fheader.adc.DC;
        case {2, 3, 4}  % Event (TTL) channels
            if ~isempty(fdata.mrk) && any(fdata.mrk(:)~=0)
                chan.markertimes = fdata.tim;
                chan.markerdata = cellstr(native2unicode(fdata.mrk(:,1)));
                if ~isvector(chan.markertimes), chan.markertimes=chan.markertimes'; end
            else
                chan.data = fdata.tim;
                if ~isvector(chan.data), chan.data=chan.data'; end
            end
        case {5}  % Marker channels
            chan.markertimes = fdata.tim;
            chan.markerdata = cellstr(native2unicode(fdata.mrk(:,1)));
        case {6}  % WaveMark channels
            error('Unsupported channel type of kind: %d', chan.kind);
        case {7}  % RealMark (Talker) channels
            if contains(lower(chan.title), 'pos')
                chan.units = 'deg';
            elseif contains(lower(chan.title), 'vel')
                chan.units = 'deg/s';
            end
            chan.samplerate = round(1/median(diff(fdata.tim)), 1);
            chan.data = single(fdata.adc);
        case {8}  % TextMark channels
            chan.markertimes = fdata.tim;
            ntimes = length(chan.markertimes);
            for jj = 1:ntimes
                chan.markerdata{end+1,1} = string(native2unicode(fdata.adc(:,jj))');
            end
        otherwise
            error('Unsupported channel type of kind: %d', chan.kind);
    end
    
    % Ensure all channel data are column vector arrays
    if isvector(chan.data) && ~iscolumn(chan.data)
        chan.data = chan.data';
    end
    
    channels{ii} = chan;
    if isempty(smr.date), smr.date=fheader.source.date; end
    smr.channelnumbers(end+1,1) = chan.number;
    smr.channelkinds(end+1,1)   = chan.kind;
    smr.channeltitles{end+1,1}  = chan.title;
end

% Ensure all channels with same samplerate are the same size
samplerates = cellfun(@(x) x.samplerate, channels);
unique_samplerates = unique(samplerates(samplerates>0.0));
if any(samplerates)
   for ii = length(unique_samplerates)
        cids = find(unique_samplerates(ii) == samplerates);
        clens = arrayfun(@(x) length(channels{x}.data), cids);
        targetlen = min(clens);
        for jj = cids'
            channels{jj}.data = channels{jj}.data(1:targetlen);
        end
    end 
end

% Set channels as a field in smr struct
smr.channels = channels;

if savefile
    % Save as a .mat file and finish import
    savefilename = [smr.filename, '.mat'];
    fullsavepath = fullfile(filepath, savefilename);
    if isfile(fullsavepath)
        dlgtxt = sprintf('%s already exists at target location.', savefilename);
        dlgtxt = [dlgtxt, '\nOverwrite existing file?'];
        dlgtitle = 'Overwrite Existing File';
        answer = questdlg(dlgtxt, dlgtitle, 'Yes', 'Cancel', 'Cancel');
        switch answer
            case 'Yes'
                % Continue and save file
            case 'Cancel'
                disp('Importing complete! File not saved.');
                return
        end
    end
    save(fullsavepath, 'smr', '-v7');
    fprintf('Importing complete! Saved as: %s\n', savefilename);
else
    disp('Importing complete! File not saved.');
end

% function dout = fillGaps(din)
% %% dout = fillGaps(din) Fill in gaps in older Spike2 recordings
% % Hannah Payne 12/16/13
% 
% dout = din;
% nsegs = size(din.data.tim,1);
% maxt = din.data.tim(end); % seconds
% samplerate = 1/prod(din.header.adc.SampleInterval);
% npoints = round(maxt*samplerate+1);
% data = zeros(1,npoints);
% tt = (1:npoints)/samplerate;
% 
% for j = 1:nsegs
%     currlength = din.header.adc.Npoints(j);
%     currdata = din.data.adc(1:currlength,j);
%     currtime = din.data.tim(j,1);
%     startind = find(tt>=currtime,1);
%     data(startind:startind+currlength-1) = currdata;
% end
% 
% dout.data.adc = data;
% dout.data.tim = [0 maxt];
% dout.header.channeltype = 'Continuous Waveform';
% end

end