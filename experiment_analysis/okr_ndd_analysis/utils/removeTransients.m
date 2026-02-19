function datout = removeTransients(datin, thresh)
% REMOVETRANSIENTS Removes transients spikes that satisfy the threshold condition

datout = datin;

% Extract left neighbors, center samples, and right neighbors
left   = datin(1:end-2);
center = datin(2:end-1);
right  = datin(3:end);

% Find spikes: large jump from both neighbors, but neighbors are close
spikes = abs(center - left)  > thresh ...
       & abs(center - right) > thresh ...
       & abs(left - right)   < thresh;

% Replace spikes with average of neighbors
datout(find(spikes) + 1) = (left(spikes) + right(spikes)) / 2;

%% TODO REMOVE ME ONCE YOU'RE CONFIDENT THE NEW VERSION ABOVE IS ROBUST
datout2 = datin;
for i = 2:length(datin)-1
    if abs(datin(i) - datin(i - 1)) > thresh && abs(datin(i) - datin(i + 1)) > thresh && abs(datin(i - 1) - datin(i + 1)) < thresh
        datout2(i) = (datin(i - 1) + datin(i + 1)) / 2;
    end
end

if ~isequal(datout, datout2)
    warning('New removeTransients method was not equal to original version!');
end


end