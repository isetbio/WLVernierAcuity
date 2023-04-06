function [P, X] = vaRGC(barOffset, params)
%% Vernier acuity experiment on cone absorptions
%
%    Testing if people can see the difference between two cases:
%      1) A straight line
%      2) Two straight lines barOffset pixels apart
%
%  Vernier acuity in human shows the positional acuity is around 6 sec of
%  arc. Here, we are analyzing stimulus, optics, and eye movements and
%  basing the calculation on absorptions.
%
% Inputs:
%   barOffset - list of bar offset
%   params    - parameter set, see s_EIParameters for more details
%
% Outputs:
%   P      - classification accuracy
%
% Note:
%   This function is adopted from s_vaAbsorptions. We removed a lot of
%   visualization code, discussion comment and print out info here.
%
% HJ/BW, ISETBIO TEAM, 2016
% JRG Added bipolar and RGC responses


%% Generate image basis
% If already computed, we load it from file. Otherwise, we make an image
% basis and save it.

% [~,imageBasis] = vaPCA(params);
nTrials = params.nTrials;

nBasis = 40;

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
    [~,alignedC] = cMosaic.compute(aligned, 'currentFlag', true, ...
        'emPaths', emPaths);
    [~,offsetC] = cMosaic.compute(offset, 'currentFlag', true, ...
        'emPaths', emPaths);    

    %%%%%%%%%%% Bipolar cell
    patchEccentricity = 2.5;
    
    cellType = {'ondiffuse','offdiffuse','onmidget','offmidget','onSBC'};
    for cellTypeInd = 1:4
        clear bpParams
        bpParams.cellType = cellType{cellTypeInd};
        
        bpParams.ecc = patchEccentricity;
        bpParams.rectifyType = 1;
        bpMosaic{cellTypeInd} = bipolar(cMosaic, bpParams);
        bpMosaic{cellTypeInd}.set('sRFcenter',1);
        bpMosaic{cellTypeInd}.set('sRFsurround',0);
        
        [~, bpNTrialsCenterTemp, bpNTrialsSurroundTemp] = bpMosaic{cellTypeInd}.compute(cMosaic,'coneTrials',alignedC);
        bpNTrialsAligned{cellTypeInd} = bpNTrialsCenterTemp-bpNTrialsSurroundTemp;
        clear bpNTrialsCenterTemp bpNTrialsSurroundTemp        
        
        [~, bpNTrialsCenterTemp, bpNTrialsSurroundTemp] = bpMosaic{cellTypeInd}.compute(cMosaic,'coneTrials',offsetC);
        bpNTrialsOffset{cellTypeInd} = bpNTrialsCenterTemp-bpNTrialsSurroundTemp;
        clear bpNTrialsCenterTemp bpNTrialsSurroundTemp
    end
    
    %%%%%%%%%% Retinal ganlion cell model
    irParams.name = 'macaque phys';
    irParams.eyeSide = 'left';
    
    % Create inner retina object
    ecc = patchEccentricity;
    irParams.eyeRadius = sqrt(sum(ecc.^2));
    irParams.eyeAngle = 0; ntrials = 0;
    innerRetina = ir(bpMosaic, irParams);
    
    mosaicParams.centerNoise = 0;
    mosaicParams.ellipseParams = [1 1 0];  % Principle, minor and theta
    % mosaicParams.axisVariance = .1;
    mosaicParams.type  = cellType;
    mosaicParams.model = 'GLM';
    
    cellType = {'on parasol','off parasol','on midget','off midget'};
    for cellTypeInd = 1:4
        mosaicParams.type = cellType{cellTypeInd};
        innerRetina.mosaicCreate(mosaicParams);
    end
    nTrials = 1; innerRetina.set('numberTrials',nTrials);
    
    % Compute the inner retina response and visualize
    disp('Computing rgc responses');
    [innerRetina, nTrialsSpikesAligned] = innerRetina.compute(bpMosaic,'bipolarTrials',bpNTrialsAligned);
    [innerRetina, nTrialsSpikesOffset]  = innerRetina.compute(bpMosaic,'bipolarTrials',bpNTrialsOffset);
       
    % trialSortedSpikesAligned = spikes2trial(nTrialsSpikesAligned);
    % trialSortedSpikesOffset = spikes2trial(nTrialsSpikesOffset);

    % Reformat the time series for the PCA analysis
    %
    % imgListX matrix contains the temporal response for a pixel in a
    % column The rows represent time samples by number of trials These are
    % the temporal responses across all trials and time points.
    
    % imgListAligned = trial2Matrix(alignedC,cMosaic);
    % imgListOffset  = trial2Matrix(offsetC,cMosaic);
    
    % imgListAligned = trial2Matrix(bpNTrialsAligned,cMosaic);
    % imgListOffset  = trial2Matrix(bpNTrialsOffset,cMosaic);
    
    imgListAligned = spikes2Matrix(nTrialsSpikesAligned);
    imgListOffset  = spikes2Matrix(nTrialsSpikesOffset);
    
    % Not-centered PCA (no demeaning, so first PC is basically the mean)
    imgList = cat(1,imgListAligned,imgListOffset);    
    
    % We make bases for each type of stimulus and then concatenate them.
    [~,~,V] = svd(imgList,'econ');
    imageBasisSpikes = V(:,1:nBasis);
    
    % Time series of weights
    weightSeries  = imgList * imageBasisSpikes;
    
    b1 = 4; b2 = 3; b3 = 2;
    figure; scatter3(weightSeries(1:194,b1),weightSeries(1:194,b2),weightSeries(1:194,b3))
    hold on; scatter3(weightSeries(194+[1:194],b1),weightSeries(194+[1:194],b2),weightSeries(194+[1:194],b3),'r')
    
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