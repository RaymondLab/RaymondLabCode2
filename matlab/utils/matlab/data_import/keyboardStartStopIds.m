function ids = keyboardStartStopIds(markers, test, train)
%KEYBOARDSTARTSTOPIDS Returns the start and end indices of experiment blocks
%   Uses regular expressions with `regeexp()` to match patterns of chars 
%   in the Keyboard channel marker data of experiment Spike2 recordings and 
%   return the start and end indices of all the matches. These regular
%   expressions usually correspond to searching for the start and end times
%   of "test" and "training" blocks within the experiment, but can also be
%   applied to experiments that do not contain training blocks at all.
%
%   If "train" is given and matches are found, the start/end indices of 
%   additional "pretest", "midtest", and "posttest" blocks will also be
%   derived from the "test" and "train" matches.
%
% Inputs:
%   markers : 1D char array of Spike2 keyboard marker data that has a
%               corresponding array of marker times.
%   test    : Regular expression to search for "test" blocks.
%   train   : (Optional) Regular expression to search for "training" blocks.
%
% Outputs:
%   ids     : Struct containing data of regexp pattern matches.
%       .test   : 2-by-N array of start/end indices of N matches for "test".
%       .ntest  : Number of matches found for regular expression "test".
%       .train  : 2-by-M array of start/end indices of M matches for "train".
%       .ntrain : Number of matches found for regular expression "train".
%       .npre   : Number of (derived) matches corresponing to "pretest" blocks.
%       .nmid   : Number of (derived) matches corresponing to "midtest" blocks.
%       .npost  : Number of (derived) matches corresponing to "posttest" blocks.
%       .all    : 2-by-(N+M) array of start/end indices of all (N+M) matches.
%       .types  : 1-by-(N+M) array of the block type of each column of ".all".
%                 Type 1 blocks are "pretests" (or "tests" if no training)
%                 Type 2 blocks are "training"
%                 Type 3 blocks are "posttests"
%                 Type 4 blocks are "midtests"
%                 Type 0 blocks correspond to an error
%
%   ----------------------------------------------------------------------
%   Author: Brian Angeles, Stanford University, 01/2025
%   ----------------------------------------------------------------------

% Ensure markers are a row char array
markers = markers(:)';

% Preallocate output struct
ids = struct;

% Get indices of the start and end of test blocks
[pA,pB] = regexp(markers(:)', test);

% Pattern is considered a fail if no test blocks are found
if isempty(pA) || isempty(pB)
    return;  % Exit early
end

ids.test = [pA; pB];
ids.ntest = length(pA);

% Get indices of start and end of all blocks
[allA,allB] = regexp(markers, [test,'|',train]);
ids.all = [allA; allB];
ids.nall = length(allA);
ids.types = zeros(ids.nall,1);

% Set all test blocks as type (1)
ids.types(ismember(allA,pA)) = 1;

if ~isempty(train)
    % Get indices (if any) of the start and end of training (2) blocks
    [tA,tB] = regexp(markers, train);
    ids.train = [tA; tB];
    ids.ntrain = length(tA);
    ids.types(ismember(allA,tA)) = 2;
else
    ids.train = [];
    ids.train = 0;
end

% Continue if there are training blocks in the recording
if ~isempty(ids.train)
    % Posttest (3) indices that occur after last training block ends
    starts = pA(pA > tB(end));
    ids.npost = length(starts);
    ids.types(ismember(allA,starts)) = 3;
    
    % Midtest (4) indices (if any) that occur during training
    starts = pA(pA < tA(end));
    ids.types(ismember(allA,starts)) = 4;

    % Pretest (1) indices that occur before first training block starts
    starts = pA(pA < tA(1));
    ids.npre = length(starts);
    ids.types(ismember(allA,starts)) = 1;
    
    % Count the remaining type 4 block types as midtest blocks
    ids.nmid = sum(ids.types == 4);
end

% Summarize the block types found (zero types are errors to investigate)
ids.uniqueTypes = unique(ids.types);

end