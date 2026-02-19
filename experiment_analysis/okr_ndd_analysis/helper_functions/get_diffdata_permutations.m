function idsout = get_diffdata_permutations(idsin)
N = length(idsin);
% Initialize pairs array
idsout = [];
% First pair: first and last element
idsout = [idsout; idsin(1), idsin(end)];
% Second pair: first and middle element (only if N > 2)
if N > 2
    midIdx = ceil(N / 2);
    idsout = [idsout; idsin(1), idsin(midIdx)];
end
% Remaining pairs: all adjacent pairs (only if N > 2)
if N > 2
    for i = 1:(N-1)
        idsout = [idsout; idsin(i), idsin(i+1)];
    end
end
% Remove duplicate pairs while preserving order
[~, uniqueIdx] = unique(idsout, 'rows', 'stable');
idsout = idsout(uniqueIdx, :);

end