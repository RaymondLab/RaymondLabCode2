function [priWC, priRad, priAmp, secWC, secRad, secAmp, ...
          pairState, frameStatus, pairAccepted, pairCost, pairDist] = ...
    trackCRPair(validProps, pairState, params)
% TRACKCRPAIR  Pair-based corneal reflection tracker using cost function.
%
% Selects both primary and secondary CRs jointly by scoring all candidate
% pairs with a weighted multi-criteria cost function. Replaces the v5
% single-CR state machine (trackCR.m).
%
% INPUTS
%   validProps  - struct array of filtered CR candidates with fields:
%                   .WeightedCentroid  [x, y]
%                   .Radius            scalar
%                   .MaxIntensity      scalar
%                   .Area              scalar
%                   .Eccentricity      scalar
%
%   pairState   - persistent state struct with fields:
%                   .prevPriPos        [x, y] — last accepted primary position
%                   .prevSecPos        [x, y] — last accepted secondary position
%                   .expectedDist      scalar — running median inter-CR distance
%                   .afterBadGap       logical — previous frame was BAD
%                   .costSlidingWin    [Nx1] — circular buffer of accepted costs
%                   .costWinCount      integer — valid entries in cost buffer
%                   .costWinIdx        integer — next write index
%                   .distSlidingWin    [Nx1] — circular buffer of accepted distances
%                   .distWinCount      integer — valid entries in distance buffer
%                   .distWinIdx        integer — next write index
%
%   params      - parameter struct with fields:
%                   .w_dist, .w_smooth, .w_geom, .w_bright  — cost weights
%                   .secConeHalfAngle  — cone half-angle (degrees)
%                   .conePenalty       — cost added for cone violation
%                   .vertCenterPenalty — cost added for vertical center violation
%                   .maxCandidates     — max shortlisted candidates
%                   .qualityGateK      — MAD multiplier for quality gate
%                   .slidingWindowSize — circular buffer capacity
%                   .minWindowCount    — warmup frames before quality gate
%                   .imgCenter         — [x, y] image center
%                   .isLeftCam         — logical
%
% OUTPUTS
%   priWC, secWC       — [x, y] primary/secondary positions
%   priRad, secRad     — radii (NaN if BAD frame)
%   priAmp, secAmp     — intensities (NaN if BAD frame)
%   pairState          — updated state struct
%   frameStatus        — uint8: 0 = GOOD, 1 = BAD
%   pairAccepted       — logical: true if real detection
%   pairCost           — total cost of selected pair (NaN if BAD)
%   pairDist           — inter-CR distance of selected pair (NaN if BAD)
%
% See also: scorePair, computeScalarWindowStats, addToScalarSlidingWindow

% -------------------------------------------------------------------------
%  Initialize outputs (default: BAD frame / fallback)
% -------------------------------------------------------------------------
priWC     = [NaN, NaN];
priRad    = NaN;
priAmp    = NaN;
secWC     = [NaN, NaN];
secRad    = NaN;
secAmp    = NaN;
frameStatus   = uint8(1);  % BAD
pairAccepted  = false;
pairCost      = NaN;
pairDist      = NaN;

nValid = numel(validProps);

% -------------------------------------------------------------------------
%  Early exit: fewer than 2 candidates — BAD frame
% -------------------------------------------------------------------------
if nValid < 2
    % Carry forward previous positions as fallback
    if ~isnan(pairState.prevPriPos(1))
        priWC = pairState.prevPriPos;
    end
    if ~isnan(pairState.prevSecPos(1))
        secWC = pairState.prevSecPos;
    end
    pairState.afterBadGap = true;
    return;
end

% -------------------------------------------------------------------------
%  Shortlist: rank by distance to prevPriPos, keep top maxCandidates
% -------------------------------------------------------------------------
dists = zeros(nValid, 1);
for c = 1:nValid
    wc = validProps(c).WeightedCentroid;
    dists(c) = sqrt((wc(1) - pairState.prevPriPos(1))^2 + ...
                    (wc(2) - pairState.prevPriPos(2))^2);
end
[~, sortOrder] = sort(dists, 'ascend');
nShort = min(nValid, params.maxCandidates);
shortlist = validProps(sortOrder(1:nShort));

% -------------------------------------------------------------------------
%  Enumerate all ordered (P, S) pairs and score each
% -------------------------------------------------------------------------
bestCost = Inf;
bestP = [];
bestS = [];
for Pi = 1:nShort
    for Si = 1:nShort
        if Pi == Si, continue; end
        cost = scorePair(shortlist(Pi), shortlist(Si), ...
            pairState.prevPriPos, pairState.prevSecPos, ...
            pairState.expectedDist, ...
            params.imgCenter, params.isLeftCam, ...
            pairState.afterBadGap, params);
        if cost < bestCost
            bestCost = cost;
            bestP = shortlist(Pi);
            bestS = shortlist(Si);
        end
    end
end

% -------------------------------------------------------------------------
%  Quality gate
% -------------------------------------------------------------------------
passQuality = false;
if pairState.costWinCount < params.minWindowCount
    % Warmup: accept unconditionally
    passQuality = true;
else
    [medCost, madCost] = computeScalarWindowStats( ...
        pairState.costSlidingWin, pairState.costWinCount, params.slidingWindowSize);
    if bestCost <= medCost + params.qualityGateK * madCost
        passQuality = true;
    end
end

% -------------------------------------------------------------------------
%  Accept or reject
% -------------------------------------------------------------------------
if passQuality
    % GOOD frame — accept this pair
    pWC = bestP.WeightedCentroid;
    sWC = bestS.WeightedCentroid;

    priWC  = pWC;
    priRad = bestP.Radius;
    priAmp = bestP.MaxIntensity;
    secWC  = sWC;
    secRad = bestS.Radius;
    secAmp = bestS.MaxIntensity;

    thisPairDist = sqrt((pWC(1) - sWC(1))^2 + (pWC(2) - sWC(2))^2);

    frameStatus  = uint8(0);  % GOOD
    pairAccepted = true;
    pairCost     = bestCost;
    pairDist     = thisPairDist;

    % Update state
    pairState.prevPriPos  = pWC;
    pairState.prevSecPos  = sWC;
    pairState.afterBadGap = false;

    % Update sliding windows
    [pairState.costSlidingWin, pairState.costWinCount, pairState.costWinIdx] = ...
        addToScalarSlidingWindow(pairState.costSlidingWin, pairState.costWinCount, ...
                                pairState.costWinIdx, bestCost, params.slidingWindowSize);
    [pairState.distSlidingWin, pairState.distWinCount, pairState.distWinIdx] = ...
        addToScalarSlidingWindow(pairState.distSlidingWin, pairState.distWinCount, ...
                                pairState.distWinIdx, thisPairDist, params.slidingWindowSize);

    % Update expected inter-CR distance from running median
    if pairState.distWinCount > 0
        [pairState.expectedDist, ~] = computeScalarWindowStats( ...
            pairState.distSlidingWin, pairState.distWinCount, params.slidingWindowSize);
    end
else
    % BAD frame — quality gate failed
    if ~isnan(pairState.prevPriPos(1))
        priWC = pairState.prevPriPos;
    end
    if ~isnan(pairState.prevSecPos(1))
        secWC = pairState.prevSecPos;
    end
    pairState.afterBadGap = true;
end

end % trackCRPair


% =========================================================================
%  LOCAL HELPER FUNCTIONS
% =========================================================================

function cost = scorePair(propP, propS, prevPriPos, prevSecPos, ...
    expectedDist, imgCenter, isLeftCam, afterBadGap, params)
% SCOREPAIR  Compute total cost for a (primary, secondary) candidate pair.

    pWC = propP.WeightedCentroid;
    sWC = propS.WeightedCentroid;

    % --- C_dist: inter-CR distance consistency ---
    pairDist = sqrt((pWC(1) - sWC(1))^2 + (pWC(2) - sWC(2))^2);
    C_dist = abs(pairDist - expectedDist) / (expectedDist + eps);

    % --- C_smooth: temporal smoothness ---
    if afterBadGap
        C_smooth = 0;
    else
        jumpPri  = sqrt((pWC(1) - prevPriPos(1))^2 + (pWC(2) - prevPriPos(2))^2);
        jumpSec  = sqrt((sWC(1) - prevSecPos(1))^2 + (sWC(2) - prevSecPos(2))^2);
        C_smooth = (jumpPri + jumpSec) / (expectedDist + eps);
    end

    % --- C_geom: geometric prior ---
    C_geom = 0;

    % Cone check: is S in the correct direction from P?
    if isLeftCam
        coneDir = [1 0];
    else
        coneDir = [-1 0];
    end
    vec = [sWC(1) - pWC(1), sWC(2) - pWC(2)];
    mag = sqrt(vec(1)^2 + vec(2)^2);
    if mag > eps
        cosAngle = dot(vec, coneDir) / mag;
        cosThresh = cosd(params.secConeHalfAngle);
        if cosAngle < cosThresh
            C_geom = C_geom + params.conePenalty;
        end
    else
        C_geom = C_geom + params.conePenalty;
    end

    % Vertical center: primary should be closer to vertical center than secondary
    priDistFromCenter = abs(pWC(2) - imgCenter(2));
    secDistFromCenter = abs(sWC(2) - imgCenter(2));
    if priDistFromCenter > secDistFromCenter
        C_geom = C_geom + params.vertCenterPenalty;
    end

    % --- C_bright: brightness ---
    C_bright = 1 / (propP.MaxIntensity + eps) + 1 / (propS.MaxIntensity + eps);

    % --- Total cost ---
    cost = params.w_dist   * C_dist   + ...
           params.w_smooth * C_smooth + ...
           params.w_geom   * C_geom   + ...
           params.w_bright * C_bright;
end


function [medVal, madVal] = computeScalarWindowStats(slidingWin, winCount, windowSize)
% COMPUTESCALARWINDOWSTATS  Compute median and MAD of a 1-D sliding window.
    n = min(winCount, windowSize);
    if n == 0
        medVal = NaN; madVal = NaN; return;
    end
    data   = slidingWin(1:n);
    medVal = median(data);
    madVal = median(abs(data - medVal));
end


function [slidingWin, winCount, winIdx] = addToScalarSlidingWindow( ...
    slidingWin, winCount, winIdx, val, windowSize)
% ADDTOSCALARSLIDINGWINDOW  Insert a scalar value into a circular sliding window.
    slidingWin(winIdx) = val;
    winIdx   = mod(winIdx, windowSize) + 1;
    winCount = min(winCount + 1, windowSize);
end
