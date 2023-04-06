function [fitProbC,fitThresh,fitParams] = PALweibullFit(levels,probC,probThresh,numTrials,fitLevels)
% PALWEIBULLFIT
%
%   [fitProbC, fitThresh, fitParams] = PALweibullFit(levels,probC,probThresh,numTrials,fitLevels)
% 
% Fits a cumulative Weibull to the (levels,probC) data. Returns the threshold at
% the probThresh value as well as the parameters needed to plot the fitted
% curve.
%
% This function requires Palamedes Toolbox 1.8
%
% Inputs:
%   stimLevels    - A N-vector representing stimuli levels corresponding to the data.
%   probC         - A N-vector of [0-1] fractions that represent performance on a task.
%   probThresh    - Probability correct called threshold.
%   numTrials     - Scalar, number of trials run for each stimulus level (assumed same for all levels).
%   fitLevels     - Stimulus levels for the fitted Weibull
%
% Returns:
%   fitProbC  - Probability correct for each of the fitLevels
%   fitThresh - Stimulus level for the probThresh
%   fitParams - Not sure how these are used, but part of PAL toolbox
%
% Example:
%    contrasts = logspace(-2,0,7); 
%    interpContrasts = logspace(-2,0,20);
%    PC = [.5 .55 0.6 0.7 0.8 0.9 1.0]; 
%    [fitPC,Thresh,paramVals] = PALweibullFit(contrasts,PC, 0.75, 200, interpContrasts);
%    vcNewGraphWin; semilogx(contrasts,PC,'o')
%    hold on; semilogx(interpContrasts,fitPC); 
%    plot(Thresh,0.75,'kx'); grid on
%
% ISETBIO Team, 2017

%% Set some parameters for the curve fitting
paramsFree     = [1 1 0 0];
numPos = round(numTrials*probC);
outOfNum       = repmat(numTrials,1,length(probC));
PF             = @PAL_Weibull;

%% Some optimization settings for the fit
options             = optimset('fminsearch');
options.TolFun      = 1e-09;
options.MaxFunEvals = 1000;
options.MaxIter     = 1000;
options.Display     = 'off';

%% Search grid specification for Palemedes
gridLevels = 100;
searchGrid.alpha = logspace(log10(levels(1)),log10(levels(end)),gridLevels);
searchGrid.beta = 10.^linspace(-2,2,gridLevels);
searchGrid.gamma = 0.5;
searchGrid.lambda = 0.0;

%% Use Palamedes grid search method
[fitParams,LL,flag] = PAL_PFML_Fit(levels(:), numPos(:), outOfNum(:), ...
            searchGrid, paramsFree, PF, 'SearchOptions', options);

%% Get threshold and deal with catastrophic cases
fitThresh = PF(fitParams, probThresh, 'inverse');
if (fitThresh < 0 || fitThresh > 1 || ~isreal(fitThresh) || isinf(fitThresh))
    fitThresh = NaN;
end

%% Provide fit psychometric function on passed stimulus levels
if (~isnan(fitThresh))
    fitProbC = PF(fitParams,fitLevels);
else
    fitProbC = NaN*ones(size(fitLevels));
end

end

