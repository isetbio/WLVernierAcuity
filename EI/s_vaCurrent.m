%% s_vaCurrent
%
%    Testing if people can see the difference between two cases:
%      1) A straight line
%      2) Two straight lines with 1 pixel apart
%
%  Vernier acuity in human shows the positional acuity is around 6 sec of
%  arc. Here, we are analyzing stimulus, optics, and eye movements and
%  basing the calculation on absorptions.
%
%  In a separate script, we will try the photocurrent.
%
% In this case we try
%    Standard retinal parameters
%    White line on a gray monitor background
%    Sweep out viewing distance
%
% HJ/BW, ISETBIO TEAM, 2016

%%
ieInit

%% Create the matched vernier stimuli.  
% Parameters are in the vaStimuli function.

[aligned, offset, scenes] = vaStimuli();
% offset.visualize;
% aligned.visualize;
% ieAddObject(offset.oiModulated); oiWindow;
% ieAddObject(scenes{2}); sceneWindow;

% Offset lines
% offsetDeg = sceneGet(scenes{1},'degrees per sample')*vparams(2).offset;
% fprintf('Offset in arc secs %.2f\n',offsetDeg*3600);

%%  Compute absorptions for multiple trials

nTrials = 100;
tSamples = aligned.length;

cMosaic = coneMosaic;

% Set the mosaic size to 15 minutes (.25 deg) because that is the spatial
% pooling size found by Westheimer and McKee
cMosaic.setSizeToFOV(0.25);

% Not sure why these have to match, but there is a bug and they do.
cMosaic.integrationTime = aligned.timeAxis(2) - aligned.timeAxis(1);  

%% For aligned or offset

tic
emPaths = zeros(nTrials, tSamples, 2);
for ii = 1 : nTrials
    cMosaic.emGenSequence(tSamples);
    emPaths(ii, :, :) = cMosaic.emPositions;
end
[alignedA, alignedC] = cMosaic.compute(aligned,...
    'emPaths',emPaths, ...
    'currentFlag',true);
toc
% cMosaic.window;

tic
emPaths = zeros(nTrials, tSamples, 2);
for ii = 1 : nTrials
    cMosaic.emGenSequence(tSamples);
    emPaths(ii, :, :) = cMosaic.emPositions;
end
[offsetA, offsetC] = cMosaic.compute(offset, ...
    'emPaths',emPaths, ...
    'currentFlag',true);
toc
% cMosaic.window;

%%  Reformat the time series for the PCA analysis

% imgListXXX matrix contains the temporal response for a pixel in a column
% The rows represent time samples by number of trials. These are the
% temporal responses across all trials and time points.

rows = cMosaic.rows;
cols = cMosaic.cols;

% Alternating between alignedC and alignedA
imgListAligned = zeros(nTrials*tSamples,rows*cols);
for tt = 1:nTrials
    lst = (1:tSamples) + (tSamples)*(tt-1);
    thisTrial = squeeze(alignedC(tt,:,:,:));
    thisTrial = permute(thisTrial,[3 1 2]);  
    thisTrial = reshape(thisTrial,tSamples,[]);
    imgListAligned(lst,:) = thisTrial;
end

imgListOffset = zeros(nTrials*tSamples,rows*cols);
for tt = 1:nTrials
    lst = (1:tSamples) + (tSamples)*(tt-1);
    thisTrial = squeeze(offsetC(tt,:,:,:));
    thisTrial = permute(thisTrial,[3 1 2]);  
    thisTrial = reshape(thisTrial,tSamples,[]);
    imgListOffset(lst,:) = thisTrial;
end

% To look at a particular trial you can set
%   cMosaic.absorptions = squeeze(offsetA(50,:,:,:));
%   cMosaic.window;
%

% Visualize the sequence
%
% imgList = abs(imgListOffset);  % Absorptions or current should be positive
% mx = max(imgListOffset(:));
% vcNewGraphWin; colormap(gray(mx));
% for ii=1:size(imgListOffset,1)
%     image(reshape(imgListOffset(ii,:),rows*cols)); 
%     drawnow; title(sprintf('%d',ii)); 
%     pause(0.05); 
% end

%% Not-centered PCA (no demeaning, so first PC is basically the mean)

% We make bases for each type of stimulus and then concatenate them.
nComponents = 5;
[~,~,V] = svd(imgListAligned,'econ');
imageBasis = V(:,1:nComponents);

[~,~,V] = svd(imgListOffset,'econ');
imageBasis = cat(2,imageBasis,V(:,1:nComponents));

% Have a look if you like
% vcNewGraphWin; colormap(gray(256))
% for ii=1:(2*nComponents)
%     imagesc(reshape(imageBasis(:,ii),rows*cols));
%     pause(0.5);
% end

imgList = cat(1,imgListAligned,imgListOffset);

% Time series of weights
weightSeries  = imgList * imageBasis;  

%% Reconstruct the input data sequence
% 
% recon = imageBasis*weights';
% vcNewGraphWin; colormap(gray(round(max(recon(:)))));
% for ii=1:size(recon,2)
%     image(reshape(recon(:,ii),rows*cols));
%     title(sprintf('%d',ii));
%     pause(0.1);
% end


%% svm classification

% TODO:  FInd a way to visualize the classifying function.
% HJ?

fprintf('SVM Classification ');

% Put the weights from each trial into the rows of a matrix
% Each row is another trial
nWeights = size(weightSeries,2);
data = zeros(2*nTrials,nWeights*tSamples);
for ii=1:(2*nTrials)
    start = (ii-1)*tSamples + 1;
    thisTrial = weightSeries(start:(start+tSamples - 1),:);
    data(ii,:) = thisTrial(:)';
end

mdl = fitcsvm(data,[ones(nTrials, 1); -ones(nTrials, 1)], 'KernelFunction', 'linear');

crossMDL = crossval(mdl);

func = @(y, yp, w, varargin) sum(abs(y(:, 1)-(yp(:, 1)>0)).*w);
classLoss = kfoldLoss(crossMDL, 'lossfun', func);
fprintf('Accuracy: %.2f%% \n', (1-classLoss) * 100);


%%
