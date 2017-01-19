function [P, mdl] = csfAbsorptions(contrasts,params)
% CSFABSORPTIONS - SVM discriminability of harmonic and uniform using absorptions
%
% The harmonic has a frequency and a contrast.  We compare that with a 0
% contrast stimulus.
%
% HJ/BW, ISETBIO TEAM, 2016
%
% See also vaAbsorptions


%% Generate image basis
% If already computed, we load it from file. Otherwise, we make an image
% basis and save it.
imageBasis = csfPCA(params);

nTrials = params.nTrials;

%% Create the aligned and offset vernier stimuli
% This could loop here on the barOffset
X = zeros(1, numel(contrasts));
P = zeros(1, numel(contrasts));
%% Compute for each contrast level?
for bb = 1:numel(contrasts)
    
    params.harmonic.contrast = contrasts(bb);
    
    fprintf('Frequency %.1f, Contrast %.3f\n',params.harmonic.freq,params.harmonic.contrast);
    [uniform, harmonic, ~, ~] = csfStimuli(params);
    % harmonic.visualize;
    
    %  Compute absorptions for multiple trials
    tSamples = uniform.length;
    
    cMosaic = coneMosaic;
    
    % Set the mosaic size to 15 minutes (.25 deg) because that is the
    % spatial pooling size found by Westheimer and McKee
    cMosaic.setSizeToFOV(params.cmFOV);
    
    % Not sure why these have to match, but there is a bug if they don't.
    cMosaic.integrationTime = uniform.timeStep;
    
    % For aligned or offset
    cMosaic.noiseFlag = 'random';
    emPaths  = cMosaic.emGenSequence(tSamples, 'nTrials', nTrials, ...
        'em', params.em);
    
    % compute absorptions for aligned and offset
    uniformA  = cMosaic.compute(uniform, 'currentFlag', false, 'emPaths', emPaths);
    harmonicA = cMosaic.compute(harmonic,'currentFlag', false, 'emPaths', emPaths);
    % cMosaic.window;
    
    % Reformat the time series for the PCA analysis
    %
    % imgListX matrix contains the temporal response for a pixel in a
    % column The rows represent time samples by number of trials These are
    % the temporal responses across all trials and time points.
    imgListUniform = trial2Matrix(uniformA,cMosaic);
    imgListharmonic  = trial2Matrix(harmonicA,cMosaic);
    
    % Not-centered PCA (no demeaning, so first PC is basically the mean)
    imgList = cat(1,imgListUniform,imgListharmonic);
    
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
    
    X(bb) = contrasts(bb);
    P(bb) = (1-classLoss) * 100;
end

end