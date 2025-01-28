function plotChannelDataUI(ax, smr, chan, varargin)
%PLOTCHANNELDATA 

% Add and validate required and optional input arguments
p = inputParser;
% validAlignMethods = {'mean', 'median', 'start', 'none'};
addRequired(p, 'ax', @mustBeNonempty);
addRequired(p, 'smr', @isstruct);
addParameter(p, 'chan', @mustBeInteger,@mustBePositive);
addParameter(p, 'xlabel', '', @mustBeText);
addParameter(p, 'ylabel', '', @mustBeText);
addParameter(p, 'xlims', [], @mustBeNumeric);
addParameter(p, 'ylims', [], @mustBeNumeric);

% Parse input arguments
parse(p, ax, smr, chan, varargin{:});

% Assign variables from parsed results
% channels    = p.Results.channels;

ydata = smr.channels{chan}.data;
xdata = single((0:length(ydata)-1)/smr.channels{chan}.samplerate);
plot(ax, xdata, ydata);
xlim(ax, [xdata(1), xdata(end)]);
end