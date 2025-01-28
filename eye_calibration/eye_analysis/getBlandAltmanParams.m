function ba = getBlandAltmanParams(avgs, diffs, saccs)
%GETBLANDALTMANPARAMS

delta_x = abs(0.3 * max(avgs));
x_min = floor(min(avgs) - delta_x);
x_max = floor(max(avgs) + delta_x);

delta_y = abs(0.35 * max(diffs));
y_min = floor(min(diffs) - delta_y);
y_max = floor(max(diffs) + delta_y);

xlims = [x_min, x_max];
ylims = [y_min, y_max];

% Compute Bias and UB/LB
b = nanmean(diffs);
UB = b + 1.96 * nanstd(diffs);
LB = b - 1.96 * nanstd(diffs);

xfit = linspace(min(avgs), max(avgs), 1000);

P1 = polyfit(avgs(~saccs), diffs(~saccs), 1);
yfit1 = polyval(P1, xfit);
txt1 = sprintf(' y = %.4fx + %.4f', P1(1), P1(2));

P3 = polyfit(avgs(~saccs), diffs(~saccs), 3);
yfit3 = polyval(P3, xfit);
txt3 = sprintf(' y = %.4fx^3 + %.4fx^2 + %.4fx + %.4f', P3(1), P3(2), P3(3), P3(4));

ba.xlims = xlims;
ba.ylims = ylims;
ba.b = b;
ba.UB = UB;
ba.LB = LB;
ba.xfit = xfit;
ba.yfit1 = yfit1;
ba.txt1 = txt1;
ba.yfit3 = yfit3;
ba.txt3 = txt3;
end