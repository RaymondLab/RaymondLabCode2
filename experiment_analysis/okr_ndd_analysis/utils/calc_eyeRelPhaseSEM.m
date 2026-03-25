function [phase,phaseSEM] = calc_eyeRelPhaseSEM(eyePhase, eyeSEM, stimPhase, stimSEM)

phase = (eyePhase - stimPhase);
phase = mod(phase, 360) - 180;
phaseSEM = sqrt((eyeSEM^2) + (stimSEM^2));

end