function [P, mdl] = vaAbsorptions(barOffset, params)
% VAABSORPTIONS - Vernier acuity experiment using cone absorptions
%
% Tests if we can use cone absorptions to detect the difference between a
% straight line vs. two straight lines offset by barOffset pixels. 
%
% The parameters for the display and stimulus are in params.vernier.
% Other control parameters are attached to params.
%
% Vernier acuity in human shows the positional acuity is around 6 sec of
% arc. Here, we are analyzing stimulus, optics, and eye movements and
% basing the calculation on absorptions.
%
% Inputs:
%   barOffset - list of bar offset
%   params    - parameter set, see s_EIParameters for more details
%
% Outputs:
%   P      - classification accuracy
%   mdl    - SVM model used to perform classification
%
% See also s_EIDefocus and many others
%
% HJ/BW, ISETBIO TEAM, 2016


%% Generate image basis
% If already computed, we load it from file. Otherwise, we make an image
% basis and save it.
imageBasis = vaPCA(params);
nTrials = params.nTrials;

%% Create the aligned and offset vernier stimuli
% This could loop here on the barOffset
X = zeros(1, numel(barOffset));
P = zeros(1, numel(barOffset));

%% Compute for each offset
for bb = 1:numel(barOffset)
    params.vernier.offset = barOffset(bb);
    fprintf('Bar offset %.1f length %.1f width %.1f fov %.2f\n',...
        barOffset(bb),params.vernier.barLength,params.vernier.barWidth,params.cmFOV);
    [aligned, offset, ~, ~] = vaStimuli(params);
    
    %  Compute absorptions for multiple trials
    tSamples = aligned.length;
    cMosaic = coneMosaic;
    
    % Sometimes we set the mosaic size to 15 minutes (.25 deg) because that is
    % the spatial pooling size found by Westheimer and McKee
    cMosaic.setSizeToFOV(params.cmFOV);
    
    % Not sure why these have to match, but there is a bug if they don't.
    cMosaic.integrationTime = aligned.timeStep;
    
    % For aligned or offset
    cMosaic.noiseFlag = 'random';
    emPaths  = cMosaic.emGenSequence(tSamples, 'nTrials', nTrials, ...
        'em', params.em);
    
    % compute absorptions for aligned and offset
    alignedA = cMosaic.compute(aligned, 'currentFlag', false, ...
        'emPaths', emPaths);
    offsetA = cMosaic.compute(offset, 'currentFlag', false, ...
        'emPaths', emPaths);
    
    % Reformat the time series for the PCA analysis
    %
    % imgListX matrix contains the temporal response for a pixel in a
    % column The rows represent time samples by number of trials These are
    % the temporal responses across all trials and time points.
    imgListAligned = trial2Matrix(alignedA,cMosaic);
    imgListOffset  = trial2Matrix(offsetA,cMosaic);
    
    % Not-centered PCA (no demeaning, so first PC is basically the mean)
    imgList = cat(1,imgListAligned,imgListOffset);
    
    % Time series of weights
    weightSeries  = imgList * imageBasis;
    
    % Start classification training
    %
    % Put the weights from each trial into the rows of a matrix
    % Each row is another trial
    nWeights = size(weightSeries,2);
    data = zeros(2*nTrials,nWeights*tSamples);
    for ii = 1 : (2*nTrials)
        start = (ii-1)*tSamples + 1;
        thisTrial = weightSeries(start:(start+tSamples - 1),:);
        data(ii,:) = thisTrial(:)';
    end
    label = [ones(nTrials, 1); -ones(nTrials, 1)];
    
    % Select some of the data (80%) as the training set.
    train_index = zeros(nTrials, 1);
    train_index(randperm(nTrials, round(0.8*nTrials))) = 1;
    train_index = train_index > 0;
    
    % The aligned and offset trials are still matched
    train_index = repmat(train_index, 2, 1);
    
    % Fit the SVM model.
    mdl = fitcsvm(data(train_index, :), label(train_index), ...
        'KernelFunction', 'linear');
    
    % predict the data not in the training set.
    yp = predict(mdl, data(~train_index, :));
    classLoss = sum(label(~train_index) ~= yp) / length(yp);
    
    X(bb) = barOffset(bb);
    P(bb) = (1-classLoss) * 100;
end

end