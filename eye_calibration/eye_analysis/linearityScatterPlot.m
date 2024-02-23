function Rsq = linearityScatterPlot(tar, vec1, vec2, keep, c)

x = vec1(keep);
y  = vec2(keep);
fit = polyfit(x,y,1);
range = linspace(min(x),max(x),100);
yfit = fit(1)*range + fit(2);
R = corrcoef(x,y);
Rsq = R(1,2).^2;

scatter(tar, vec2(keep), vec1(keep), 4, c(keep), '.');
plot(tar, yfit, range, 'k', 'lineWidth',2);
colormap(tar, hsv);

end