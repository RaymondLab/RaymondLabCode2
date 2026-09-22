function [r2, maxr2, maxr2_location] = linearityAlign(x, y)
%LINEARITYALIGN

% If input vectors are not the same length, shorten the longer one
if length(x) < length(y)
    y = y(1:length(x));
elseif length(x) > length(y)
    x = x(1:length(y));    
end

r2 = nan(3901,1);

for i = -4000:-100              
    if i < 0
        x2 = x((-i)+1:end);
        y2 = y(1:end+i);
    elseif i == 0
        x2 = x;
        y2 = y;
    else
        x2 = x(1:end-i);
        y2 = y((i)+1:end);
    end
    R2 = corrcoef(x2,y2);
    r2(i+4001) = R2(1,2).^2;
end

maxr2 = max(r2(:,1));
maxr2_location = find(r2(:,1) == max(r2(:,1)),1) - 4000;