function [goodStarts,goodStops,returnedError] = getSaccadeStartEndTimes(saccades)
%GETSACCADESTARTENDTIMES

returnedError = 0;

% Saccade start and end times
sacEndPoints = diff(~saccades);
goodStarts = find(sacEndPoints == 1);
goodStops = find(sacEndPoints == -1);

if goodStarts(1) > goodStops(1)
    goodStarts = [1; goodStarts];
elseif goodStarts(end) > goodStops(end)
    goodStops = [goodStops; length(~saccades)];
elseif isempty(goodStarts) && isempty(goodStops)
    returnedError = 1;
    disp('No saccades: Cannot calculate piecewise linearity')
    return
end

end