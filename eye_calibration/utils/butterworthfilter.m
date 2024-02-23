function [dataFiltered] = butterworthfilter(data, cutoff, samplerate, N)
%BUTTERWORTHFILTER Applies N-th order low-pass Butterworth filter on data.

if ~exist('N', 'var')
    N = 9;
end

[z, p, k] = butter(N, cutoff/(samplerate/2), 'low');
[sos, g] = zp2sos(z, p, k);

dataFiltered = filtfilt(sos, g, data);

end