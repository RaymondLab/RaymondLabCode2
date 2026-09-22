function B = butterworthfilter(A, fc, samplerate, varargin)
%BUTTERWORTHFILTER Butterworth filter
%   Applies an n-th order, zero-phase, `ftype` Butterworth filter with 
%   cutoff frequency `fc` on given array `A`.
%
%   If fc is a two-element vector, fc = [fc1 fc2], BUTTERWORTHFILTER 
%   returns an order 2n bandpass filter with passband  fc1 < fc < fc2.
%   B = BUTTERWORTHFILTER(A,fc,'high') designs a highpass filter.
%   B = BUTTERWORTHFILTER(A,fc,'low') designs a lowpass filter.
%   B = BUTTERWORTHFILTER(A,fc,'stop') is a bandstop filter if fc = [fc1 fc2].
%
% Inputs:
%   A          : 1D array of samples to apply filter
%   fc         : cutoff frequency for low-pass filter
%   samplerate : samplerate corresponding to array `a` (default 1000)
%   n          : filter order (default 9)
%   ftype      : filter type (default "low")
%   
% Outputs:
%   B          : low-pass filtered array
%
%   ----------------------------------------------------------------------
%   Author: Brian Angeles, Stanford University, 03/2025
%   ----------------------------------------------------------------------

% Add and validate required and optional input arguments
p = inputParser;
ftypes = {'low', 'high', 'stop'};
addRequired(p, 'A', @(x)~isempty(x)&&isnumeric(x));
addRequired(p, 'fc', @(x)~isempty(x));
addRequired(p, 'samplerate', @(x)isnumeric(x)&&(x>0));
addParameter(p, 'n', 9, @(x)isnumeric(x)&&(x>0));
addParameter(p, 'ftype', "low", @(x)any(validatestring(x,ftypes)));

% Parse input arguments
parse(p, A, fc, samplerate, varargin{:});
n = p.Results.n;
ftype = p.Results.ftype;

nq = samplerate / 2;  % Half of sample rate (nyquist)
[z,p,k] = butter(n, fc/nq, ftype);
[sos, g] = zp2sos(z, p, k);
B = filtfilt(sos, g, A);
end