function [bestFit,yFit,bestGof,dispName] = fitBestModel(x, y, varargin)
% FITBESTMODEL Fits multiple models to time series data and selects the best.
%
%   [bestFit, yFit, bestGof, dispName] = fitBestModel(x, y)
%   [bestFit, yFit, bestGof, dispName] = fitBestModel(x, y, 'Verbose', true)
%
%   Inputs:
%       x           - Independent variable (e.g., block number), vector
%       y           - Dependent variable (e.g., metric values), vector
%
%   Optional Name-Value Pairs:
%       'Verbose'              - Print comparison table to command window
%                                (default: true)
%       'Models'               - Cell array of fit types to test
%                                (default: {'poly1','poly2','exp1','exp2','power1','log'})
%       'MinRSquared'          - R² threshold below which the fit is flagged
%                                as not meaningful in dispName (default: 0.1)
%       'PreferSimpleWhenTied' - If true, prefer 'poly1' whenever its AIC is
%                                within AICTolerance of the lowest-AIC model
%                                (default: true)
%       'AICTolerance'         - ΔAIC threshold for the parsimony rule above.
%                                Per Burnham & Anderson, ΔAIC < AICTolerance indicates
%                                no meaningful evidence against the simpler
%                                model (default: 3.5)
%
%   Outputs:
%       bestFit     - cfit object for the best model
%       yFit        - Fitted y-values evaluated at input x (same size as x)
%       bestGof     - goodness-of-fit struct for the best model
%       dispName    - Sprintf string of the equation, ready for DisplayName.
%                     If R² < MinRSquared, this is replaced with a
%                     "no meaningful trend" message instead of the equation.
%
%   Example:
%       x = (1:60)';
%       y = 3*exp(-0.05*x) + 0.2*randn(60,1);
%       [bestFit, yFit, bestGof, dispName] = fitBestModel(x, y);
%       plot(x, y, 'o'); hold on;
%       plot(x, yFit, '-', 'DisplayName', dispName);
%       legend('show');
%
%   Notes:
%       - Model selection uses AIC (Akaike Information Criterion), with an
%         optional parsimony tiebreaker that defaults to 'poly1' when the
%         AIC gap to the nominal best model is below AICTolerance.
%       - Fits use Bisquare robust weighting where supported, so that
%         single leverage points cannot dominate model selection.
%       - Models that fail to converge are skipped and noted in the output.
%       - Requires the Curve Fitting Toolbox.

% --- Parse inputs ---
p = inputParser;
addRequired(p, 'x', @isnumeric);
addRequired(p, 'y', @isnumeric);
addParameter(p, 'Verbose', true, @islogical);
addParameter(p, 'Models', {'poly1','poly2','exp1','exp2','power1','log'}, @iscell);
addParameter(p, 'MinRSquared', 0.1, @(v) isnumeric(v) && isscalar(v));
addParameter(p, 'PreferSimpleWhenTied', true, @islogical);
addParameter(p, 'AICTolerance', 5, @(v) isnumeric(v) && isscalar(v) && v >= 0);
parse(p, x, y, varargin{:});

verbose              = p.Results.Verbose;
modelList            = p.Results.Models;
minRSquared          = p.Results.MinRSquared;
preferSimpleWhenTied = p.Results.PreferSimpleWhenTied;
aicTolerance         = p.Results.AICTolerance;

% --- Store original x for output evaluation ---
xOrig = x(:);
inputWasRow = isrow(x);

% --- Ensure column vectors ---
x = x(:);
y = y(:);

% --- Validate inputs ---
if length(x) ~= length(y)
    error('fitBestModel:dimMismatch', 'x and y must have the same length.');
end

% Remove NaN/Inf pairs
validIdx = isfinite(x) & isfinite(y);
if sum(validIdx) < length(x)
    warning('fitBestModel:removedNaN', ...
        'Removed %d non-finite data points.', sum(~validIdx));
end
xClean = x(validIdx);
yClean = y(validIdx);
n = length(xClean);

if n < 3
    error('fitBestModel:tooFewPoints', 'Need at least 3 valid data points.');
end

% --- Storage ---
nModels    = length(modelList);
fitResults = cell(nModels, 1);
gofResults = cell(nModels, 1);
aicValues  = nan(nModels, 1);
fitSuccess = false(nModels, 1);

% --- Number of coefficients per built-in fit type ---
numCoeffsMap = containers.Map( ...
    {'poly1','poly2','poly3','poly4','poly5','poly6','poly7','poly8','poly9', ...
     'exp1','exp2','power1','power2','log'}, ...
    {2, 3, 4, 5, 6, 7, 8, 9, 10, ...
     2, 4, 2, 4, 2});

% --- Fit each model ---
for i = 1:nModels
    modelName = modelList{i};
    try
        opts = fitoptions(modelName);
        % MaxIter / MaxFunEvals only exist for NonlinearLeastSquares; linear
        % fits (poly*, log) use LinearLeastSquares and would error on them.
        optProps = properties(opts);
        if any(strcmp(optProps, 'MaxIter'));     opts.MaxIter = 1000;     end
        if any(strcmp(optProps, 'MaxFunEvals')); opts.MaxFunEvals = 2000; end
        % Down-weight outliers so a single leverage point cannot drive selection.
        if any(strcmp(optProps, 'Robust'));      opts.Robust = 'Bisquare'; end

        [fitResults{i}, gofResults{i}] = fit(xClean, yClean, modelName, opts);

        if isKey(numCoeffsMap, modelName)
            k = numCoeffsMap(modelName);
        else
            k = numcoeffs(fitResults{i});
        end

        sse = gofResults{i}.sse;
        if sse > 0 && n > k
            aicValues(i) = n * log(sse / n) + 2 * k;
        else
            aicValues(i) = Inf;
        end

        fitSuccess(i) = true;

    catch ME
        if verbose
            fprintf('  %-10s  FAILED: %s\n', modelName, ME.message);
        end
    end
end

% --- Select best model by AIC ---
if ~any(fitSuccess)
    error('fitBestModel:allFailed', 'All model fits failed.');
end

[minAIC, bestIdx] = min(aicValues);

% --- Parsimony tiebreaker: prefer poly1 when within AICTolerance ---
% AIC differences below ~2 are not considered meaningful evidence against
% the simpler model (Burnham & Anderson). When that's the case, fall back
% to the most interpretable model (poly1) instead of letting tiny SSE
% differences pick a curved fit.
parsimonyApplied = false;
if preferSimpleWhenTied
    poly1Idx = find(strcmp(modelList, 'poly1'), 1);
    if ~isempty(poly1Idx) && fitSuccess(poly1Idx) && poly1Idx ~= bestIdx
        if aicValues(poly1Idx) - minAIC <= aicTolerance
            bestIdx = poly1Idx;
            parsimonyApplied = true;
        end
    end
end

bestFit = fitResults{bestIdx};
bestGof = gofResults{bestIdx};
bestName = modelList{bestIdx};

% --- Evaluate fit at original x values ---
yFit = feval(bestFit, xOrig);
if inputWasRow
    yFit = yFit';
end

% --- Build DisplayName string ---
if bestGof.rsquare < minRSquared
    dispName = sprintf('no meaningful trend (R^2=%.3f)', bestGof.rsquare);
else
    dispName = getDisplayName(bestFit, bestName, bestGof);
end

% --- Verbose: print comparison table ---
if verbose
    fprintf('\n  %-10s  %8s  %8s  %10s  %8s\n', ...
        'Model', 'R²', 'Adj R²', 'AIC', 'RMSE');
    fprintf('  %s\n', repmat('-', 1, 52));
    for i = 1:nModels
        if fitSuccess(i)
            tag = '';
            if i == bestIdx
                if parsimonyApplied
                    tag = '  <-- selected (parsimony)';
                else
                    tag = '  <-- best';
                end
            elseif aicValues(i) == minAIC
                tag = '  <-- min AIC';
            end
            fprintf('  %-10s  %8.4f  %8.4f  %10.2f  %8.4f%s\n', ...
                modelList{i}, ...
                gofResults{i}.rsquare, gofResults{i}.adjrsquare, ...
                aicValues(i), gofResults{i}.rmse, tag);
        end
    end
    fprintf('  %s\n', repmat('-', 1, 52));
    if parsimonyApplied
        fprintf('  poly1 selected over min-AIC model (ΔAIC = %.2f <= %.2f).\n', ...
            aicValues(bestIdx) - minAIC, aicTolerance);
    end
    fprintf('  DisplayName: %s\n\n', dispName);
end

end

%% ========================================================================
%  LOCAL FUNCTION: Build DisplayName string (equation + R²)
%  ========================================================================
function dispName = getDisplayName(f, modelName, gof)
    coeffs = coeffvalues(f);

    switch modelName
        case 'poly1'
            eqn = sprintf('y = %.3gx + %.3g', coeffs(1), coeffs(2));
        case 'poly2'
            eqn = sprintf('y = %.3gx^2 + %.3gx + %.3g', ...
                coeffs(1), coeffs(2), coeffs(3));
        case 'exp1'
            eqn = sprintf('y = %.3ge^{%.3gx}', coeffs(1), coeffs(2));
        case 'exp2'
            eqn = sprintf('y = %.3ge^{%.3gx} + %.3ge^{%.3gx}', ...
                coeffs(1), coeffs(2), coeffs(3), coeffs(4));
        case 'power1'
            eqn = sprintf('y = %.3gx^{%.3g}', coeffs(1), coeffs(2));
        case 'power2'
            eqn = sprintf('y = %.3gx^{%.3g} + %.3g', ...
                coeffs(1), coeffs(2), coeffs(3));
        case 'log'
            eqn = sprintf('y = %.3gln(x) + %.3g', coeffs(1), coeffs(2));
        otherwise
            cnames = coeffnames(f);
            parts = cell(1, length(cnames));
            for j = 1:length(cnames)
                parts{j} = sprintf('%s=%.3g', cnames{j}, coeffs(j));
            end
            eqn = strjoin(parts, ', ');
    end

    dispName = sprintf('%s (R^2=%.3f)', eqn, gof.rsquare);
end