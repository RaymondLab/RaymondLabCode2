function vidAmp = vidTimeFcn_APP(app, tvid, vidH, freq, tscale)
    varsTemp = [sin(2*pi*freq*tvid(:)/tscale) cos(2*pi*freq*tvid(:)/tscale) ones(size(tvid(:)))];
    yTemp = [0; diff(vidH)]/mean(diff(tvid));
    mask = abs(yTemp)<100;
    bVidTemp = regress(yTemp(mask), varsTemp(mask,:));
    vidAmp = -sqrt(bVidTemp(1)^2+bVidTemp(2)^2);
end