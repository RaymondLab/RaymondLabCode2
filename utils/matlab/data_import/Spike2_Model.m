classdef Spike2_Model < handle
    %SPIKE2_MODEL Data model for Spike2 recording data

    % Define public events that can be used to trigger outside functions
    events
        dataChanged
        filterValueChanged
    end

    % Observed public properties
    properties (SetObservable)
        % Spike2 struct file generated via `importSpike2File` function
        smr (1,1) struct
        % Sample rate of derived raw eye position signal in Hz
        samplerate (1,1) double {mustBeNumeric,mustBePositive} = 1000.0;
        % Frequency of stimulus
        stimulusFrequency (1,1) double {mustBeNumeric,mustBePositive} = 1.0;
        % Eye calibration scaling factors of respective magnet channels
        scaleCh (2,1) double {mustBeNumeric} = [NaN; NaN];
        % Cutoff frequency (Hz) used for lowpass Butterworth filter
        lowpassCutoff (1,1) double {mustBeNumeric,mustBePositive} = 40.0;  %Todo: 15?
        % Window (ms) used for movingslope ("savgol") filter
        savgolWindow (1,1) double {mustBeInteger,mustBePositive} = 30;  %Todo: 11?
    end

    % Private properties
    properties (SetAccess = protected)
        % Struct containing the indices of channels in smr property
        chid (1,1) struct
        % Struct containing the samplerates of key signal channels
        fs (1,1) struct
        % Timetable of raw and processed recording traces
        data (:,:) timetable
        % Keyboard marker channel data
        kbdata (:,1) char
        % Keyboard channel marker times
        kbtimes (:,1) double {mustBeNumeric}
        % TextMark marker channel data
        tmdata (:,1) string
        % TextMark channel marker times
        tmtimes (:,1) double {mustBeNumeric}
        % Indices of the start/end of each block of experiment
        blockids (2,:) double {mustBeNumeric}
        % Indices of the start/end of each block of experiment
        blocktimes (:,2) double {mustBeNumeric}
        % Corresponding type of each block
        blocktypes (:,1) double {mustBeNumeric}
        % Total number of blocks in experiment
        nblocks (1,1) double {mustBeInteger}
        % Index of schema successfully used for block search
        schemaidx (1,1) double {mustBeInteger}
        % Flag to bypass listeners when updating observed properties
        isUpdating (1,1) logical = false;
    end

    methods
        function obj = Spike2_Model()
            %SPIKE2_MODEL Class constructor
            addlistener(obj, 'smr', 'PostSet',@obj.change_smr);
            addlistener(obj, 'scaleCh', 'PostSet',@obj.change_scaleCh);
            addlistener(obj, 'lowpassCutoff', 'PostSet',@obj.change_lowpassCutoff);
            addlistener(obj, 'savgolWindow', 'PostSet',@obj.change_savgolWindow);
        end
        
        % Retrieve and set the start/end indices of blocks via keyboard data
        function setBlockIds(obj, test, train)
            if ~exist('test','var'), test=''; end
            if ~exist('train','var'), train=''; end
            blocks = blockSchemaSearch(obj.kbdata, test, train);
            obj.schemaidx = blocks.schemaIdx;
            obj.nblocks = blocks.nall;
            obj.blockids = blocks.all;
            obj.blocktimes = obj.kbtimes(obj.blockids)';
            obj.blocktypes = blocks.types;
        end
        
        % Define setter functions to set properties from outside the class
        function setSmr(obj, data)
            obj.smr = data;
        end

        function setStimulusFrequency(obj, value)
            obj.stimulusFrequency = value;
        end

        function setScaleCh(obj, value)
            if isrow(value)
                value = value';
            end
            obj.scaleCh = value;
        end

        function setLowpassCutoff(obj, value)
            obj.lowpassCutoff = value;
        end

        function setSavgolWindow(obj, value)
            obj.savgolWindow = round(value);
        end
    end

    methods (Access = protected)
        function change_lowpassCutoff(obj, varargin)
            if isempty(obj.data) || obj.isUpdating, return; end
            obj.isUpdating = true;
            obj.data.hepos = butterworthfilter(obj.data.heposraw, obj.lowpassCutoff, obj.samplerate);
            obj.data.hevel = movingslope(obj.data.hepos, obj.savgolWindow) * obj.samplerate;
            % Notify listeners that the lowpass cutoff value was changed
            notify(obj, 'filterValueChanged');
            obj.isUpdating = false;
        end

        function change_savgolWindow(obj, varargin)
            if isempty(obj.data) || obj.isUpdating, return; end
            obj.isUpdating = true;
            obj.data.hevel = movingslope(obj.data.hepos, obj.savgolWindow) * obj.samplerate;
            % Notify listeners that the savgol window value was changed
            notify(obj, 'filterValueChanged');
            obj.isUpdating = false;
        end

        function change_scaleCh(obj, varargin)
            if isempty(obj.smr) || obj.isUpdating, return; end
            obj.isUpdating = true;
            hepos1 = obj.smr.channels{obj.chid.hepos1}.data * obj.scaleCh(1);
            hepos2 = obj.smr.channels{obj.chid.hepos2}.data * obj.scaleCh(2);
            tt = timetable(hepos1+hepos2, ...
                'SampleRate',obj.samplerate, ...
                'VariableNames',{'heposraw'});
            tt.hevelraw = movingslope(tt.heposraw, 3) * obj.samplerate;
            tt.heposdefault = butterworthfilter(tt.heposraw, 40.0, obj.samplerate); 
            tt.heveldefault = movingslope(tt.heposdefault, 3) * obj.samplerate;
            tt.hepos = butterworthfilter(tt.heposraw, obj.lowpassCutoff, obj.samplerate);
            tt.hevel = movingslope(tt.hepos, obj.savgolWindow) * obj.samplerate;
            try tt.htvel = obj.smr.channels{obj.chid.htvel}.data; end
            try tt.hhvel = obj.smr.channels{obj.chid.hhvel}.data; end
            try tt.HTVEL = obj.smr.channels{obj.chid.HTVEL}.data; end
            try tt.HHVEL = obj.smr.channels{obj.chid.HHVEL}.data; end
            obj.data = tt;
            obj.isUpdating = false;
            % Notify listeners that the data property was changed
            notify(obj, 'dataChanged');
        end

        function change_smr(obj, varargin)
            if obj.isUpdating, return; end
            obj.isUpdating = true;
            % Redefine the cid struct based on the new smr file
            channeltitles = obj.smr.channeltitles;
            obj.chid = struct;
            obj.chid.hepos1 = find(strcmpi(channeltitles, 'hepos1'));
            obj.chid.hepos2 = find(strcmpi(channeltitles, 'hepos2'));
            obj.chid.htvel = find(ismember(channeltitles, 'htvel'));
            obj.chid.hhvel = find(ismember(channeltitles, 'hhvel'));
            obj.chid.HTVEL = find(ismember(channeltitles, 'HTVEL'));
            obj.chid.HHVEL = find(ismember(channeltitles, 'HHVEL'));
            obj.chid.keyboard = find(strcmpi(channeltitles, 'keyboard'));
            obj.chid.textmark = find(strcmpi(channeltitles, 'textmark'));
            % Update Keyboard data and times
            if isscalar(obj.chid.keyboard)
                kbidx = obj.chid.keyboard;
                kbmd = cell2mat(obj.smr.channels{kbidx}.markerdata);
                obj.kbdata = kbmd(:);
                kbmt = obj.smr.channels{kbidx}.markertimes;
                obj.kbtimes = kbmt(:);
            end
            % Update (if available) TextMark data and times
            if isscalar(obj.chid.textmark)
                tmidx = obj.chid.textmark;
                tmmd = obj.smr.channels{tmidx}.markerdata;
                tmmd = cellfun(@(x) regexprep(x,'\o{0}',''), tmmd);
                obj.tmdata = tmmd(:);
                tmmt = obj.smr.channels{tmidx}.markertimes;
                obj.tmtimes = tmmt(:);
            end
            % Try detecting the start/end of blocks in the recording
            obj.setBlockIds();
            % Reset the scaling factors to NaN values
            obj.scaleCh = [NaN; NaN];
            % Update samplerate property
            obj.fs.hepos1 = obj.smr.channels{obj.chid.hepos1}.samplerate;
            obj.fs.hepos2 = obj.smr.channels{obj.chid.hepos2}.samplerate;
            % Todo: what do we do about slightly different samplerates?
            if obj.fs.hepos1 ~= obj.fs.hepos2 
                obj.samplerate = double((obj.fs.hepos1 + obj.fs.hepos1) / 2.0);
            else
                obj.samplerate = obj.fs.hepos1;
            end
            obj.isUpdating = false;
        end
    end
    
end