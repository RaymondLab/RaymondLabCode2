function [cyclemat,cyclemean] = getCycleData(a, cycleids, cycleperiod)
%GETCYLEDATA Summary of this function goes here

ncycles = length(cycleids)-1;
cyclemat = nan(cycleperiod, ncycles);
cyclemean = nan;
for ii = 1:ncycles
    start = cycleids(ii);
    stop = cycleids(ii+1)-1;
    acycle = a(start:stop);
    minperiod = min([cycleperiod, length(acycle)]);
    cyclemat(1:minperiod,ii) = acycle(1:minperiod);
end
cyclemean = nanmean(cyclemat, 2);

end

