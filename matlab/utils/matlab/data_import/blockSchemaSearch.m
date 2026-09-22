function blocks = blockSchemaSearch(markers, test, train)
%BLOCKSCHEMASEARCH Attempts various schemas to return block start/end indices
%   Detailed explanation goes here
%
% Inputs:
%   markers : 1D char array of Spike2 keyboard marker data that has a
%               corresponding array of marker times.
%   test    : (Optional) Regular expression to search for "test" blocks.
%   train   : (Optional) Regular expression to search for "training" blocks.
%
%   ----------------------------------------------------------------------
%   Author: Brian Angeles, Stanford University, 01/2025
%   ----------------------------------------------------------------------

% Define the pattern schemas to try
TestTrainSchemas = { ...
    {'P[^p]*p'; 'T[^t]*t'}; ...                             % Newest protocols
    {'1XS[^lIs]*lIs|2XS[^3x]*3x'; '3xXL[^0xl]*0xl'}; ...    % Trace Consolidation & Amin VOR
    {'4X[^1x]*1x|4X[^0x]*0x'; '3X[^1x]*1x|3X[^0x]*0x'}; ... % Amin OKR
    {'1S[^0p]*0p'; 'LS[^0p]*0p'}                            % Sriram OKR
};

% Ensure markers are a row char array
markers = markers(:)';

blocks = struct;
if ~isempty(test)
    % Use custom test/train schema if arguments are provided
    blocks = keyboardStartStopIds(markers, test, train);
    blocks.schemaIdx = 0;
else
    % Otherwise, try each schema and return the first successful attempt
    for ii = 1:length(TestTrainSchemas)
        testSchema = TestTrainSchemas{ii}{1};
        trainSchema = TestTrainSchemas{ii}{2};
        blocks = keyboardStartStopIds(markers, testSchema, trainSchema);
        if ~isempty(fieldnames(blocks)) & ~isempty(blocks.test)
            fprintf('%d blocks found!\n', blocks.nall);
            blocks.schemaIdx = ii;
            return % Exit early on successful attempt
        else
            fprintf('No blocks found. Trying next schema...\n');
        end
    end
end