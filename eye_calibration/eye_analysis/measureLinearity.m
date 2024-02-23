function d = measureLinearity(vec1, vec2, keep)

if ~exist('keep', 'var')
    keep = ones(length(vec1),1);
end

d.fit = polyfit(vec1,vec2,1);

d.range = linspace(min(vec1),max(vec1),100);
d.range_Norm = d.range - d.range(1);

d.yfit = d.fit(1)*d.range + d.fit(2);
d.yfit_Norm = d.yfit - d.yfit(1);

d.R = corrcoef(vec1,vec2);
d.Rsq = d.R(1,2).^2;
d.slope = d.fit(1);


