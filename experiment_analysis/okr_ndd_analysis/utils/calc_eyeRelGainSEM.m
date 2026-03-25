function [gain,gainSEM] = calc_eyeRelGainSEM(eyeAmp, eyeSEM, stimAmp, stimSEM)
    
gain = eyeAmp / stimAmp;
a_sq = (eyeSEM / eyeAmp)^2;
b_sq = (stimSEM / stimAmp)^2;
gainSEM = gain * sqrt(a_sq + b_sq);

end