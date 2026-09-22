function [r2, maxr2, maxr2_location] = linearityAlignV2(x, y)
%LINEARITYALIGNV2
[r2, lags] = xcorr(y, x);
ccnorm = r2 / max(r2);
[~, maxr2idx] = max(ccnorm);
r2 = r2.^2;
maxr2 = r2(maxr2idx);
maxr2_location = lags(maxr2idx);
end