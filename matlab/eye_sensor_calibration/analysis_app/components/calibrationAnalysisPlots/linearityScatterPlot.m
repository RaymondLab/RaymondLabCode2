function Rsq = linearityScatterPlot(tar, vec1, vec2, keep, c)

x = vec1(keep);
y = vec2(keep);
llimits = min(min(x), min(y));
ulimits = max(max(x), max(y));
fit = polyfit(x, y, 1);
range = linspace(llimits, ulimits, 100);
yfit = fit(1)*range + fit(2);
R = corrcoef(x, y);
Rsq = R(1, 2).^2;

scatter(tar, vec1(keep), vec2(keep), 4, c(keep), '.');
plot(tar, range, yfit, 'k', 'lineWidth',2);
colormap(tar, hsv);
ylim(tar, [llimits,ulimits]);
xlim(tar, [llimits,ulimits]);

end