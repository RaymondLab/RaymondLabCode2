function [amplitude, phase, freq, fitTrace, r2] = fit_sineWave(data, samplerate, freq)
%FIT_SINEWAVE

% set up thing to fit
segLength = length(data);
segTime = (1:segLength)/samplerate; 

% if no frequency given, estimate
if ~exist('freq', 'var')
    
    Y = fft(data(~isnan(data)));
    P2 = abs(Y/segLength);
    powerTrace = P2(1:segLength/2+1);
    powerTrace(2:end-1) = 2*powerTrace(2:end-1);
    freqDomain = samplerate*(0:(segLength/2))/segLength;
    
    minFreq = .1;
    powerTrace(freqDomain < minFreq) = [];
    freqDomain(freqDomain < minFreq) = [];
    
    pks = findpeaks(powerTrace, samplerate);
    biggestPeak = max(pks);
    exactFreq = freqDomain(find(powerTrace == biggestPeak));
    roundedFreq = round(exactFreq, 1);
    freq = roundedFreq;
end

% generate vars
y1 = sin(2*pi*freq*segTime(:));
y2 = cos(2*pi*freq*segTime(:));
constant = ones(segLength,1);
vars = [y1 y2 constant];

% Do regression, calculate amplitude, phase, r^2, and fit
[b,~,~,~,stat] = regress(data, vars);
amplitude = sqrt(b(1)^2+b(2)^2);
phase = rad2deg(atan2(b(2), b(1)));
r2 = stat(1);
fitTrace = sin(2*pi*freq*segTime + deg2rad(phase))*amplitude;

end