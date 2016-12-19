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

nTrials = 250;
tStep   = 5;   % ms

% Set parameters for the vernier stimulus
clear p; 
p.barOffset = 3;     % Pixels on the display
p.barWidth  = 3;     % Pixels on the display
p.tsamples = (-60:tStep:70)*1e-3;   % In second
p.timesd = 20*1e-3;                 % In seconds

%% Create the matched vernier stimuli

[aligned, offset, scenes,tseries] = vaStimuli(p);
% offset.visualize;
% aligned.visualize;
% ieAddObject(offset.oiModulated); oiWindow;
% ieAddObject(scenes{2}); sceneWindow;
% vcNewGraphWin; plot(p.tsamples,tseries)

% Offset lines
% offsetDeg = sceneGet(scenes{1},'degrees per sample')*vparams(2).offset;
% fprintf('Offset in arc secs %.2f\n',offsetDeg*3600);

%%  Compute absorptions for multiple trials

tSamples = aligned.length;
cMosaic = coneMosaic;

% Set the mosaic size to 15 minutes (.25 deg) because that is the spatial
% pooling size found by Westheimer and McKee
cMosaic.setSizeToFOV(0.25);

% Not sure why these have to match, but there is a bug and they do.
cMosaic.integrationTime = aligned.timeStep;  

%% For aligned or offset
cMosaic.os.noiseFlag = 'random';
cMosaic.noiseFlag    = 'random';

emPaths = cMosaic.emGenSequence(tSamples,'nTrials',nTrials);

tic
[~, alignedC] = cMosaic.compute(aligned, ...
    'emPaths',emPaths, ...
    'currentFlag',true);
toc
% cMosaic.current = squeeze(mean(alignedC,1));
% cMosaic.window;

tic
[~, offsetC] = cMosaic.compute(offset, ...
    'emPaths',emPaths, ...
    'currentFlag',true);
toc
% cMosaic.current = squeeze(mean(offsetC,1));
% cMosaic.window;

%%  Reformat the time series for the PCA analysis

% imgListXXX matrix contains the temporal response for a pixel in a column
% The rows represent time samples by number of trials. These are the
% temporal responses across all trials and time points.

imgListAligned = trial2Matrix(alignedC,cMosaic);
imgListOffset  = trial2Matrix(offsetC,cMosaic);

%% Not-centered PCA (no demeaning, so first PC is basically the mean)

% We make bases for each type of stimulus and then concatenate them.
nComponents = 5;
[~,~,V] = svd(imgListAligned,'econ');
imageBasis = V(:,1:nComponents);

[~,~,V] = svd(imgListOffset,'econ');
imageBasis = cat(2,imageBasis,V(:,1:nComponents));

% imgList = cat(1,imgListAligned,imgListOffset);
imgList = cat(1,imgListAligned,imgListAligned);

% Time series of weights
weightSeries  = imgList * imageBasis;


%% svm classification

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

%
% func = @(y, yp, w, varargin) sum(abs(y(:, 1)-(yp(:, 1)>0)).*w);
% classLoss = kfoldLoss(crossMDL, 'lossfun', func);

% train with 80% of data
label = [ones(nTrials, 1); -ones(nTrials, 1)];
train_index = zeros(2*nTrials, 1);
train_index(randperm(2*nTrials, 0.8*2*nTrials)) = 1;
train_index = train_index > 0;

mdl = fitcsvm(data(train_index, :), label(train_index), ...
    'KernelFunction', 'linear');

% predict on test set
yp = predict(mdl, data(~train_index, :));
classLoss = sum(label(~train_index) ~= yp) / length(yp);
fprintf('Accuracy: %.2f%% \n', (1-classLoss) * 100);

%% Visualize the classification function

% Field of view of the mosaic
% Eye movement pattern
% Stimulus offset, more ...

% See 'Support Vector Machines for Binary Classification' in fitcsvm
%
%  mdl.Beta appears to be the term that we multiply as an inner product
%  with the data to determine whether we are in type A or not A.
%
%  x'*beta + Bias, I think.
%
%  When we get mdl.Beta, it has size of (2*nComponents)*tSamples.   So, if
%  there are 150 tSamples and 10 spatial image components, then Beta is
%  1500.
%
%  I think that on a single trial we have 150*10 numbers, which we derive
%  by taking the inner product of the 10 spatial basis functions on every
%  temporal sample.  If we want to express the 1500 numbers as a time
%  series, we would multiple the 10 weights times the spatial image basis
%  at each of the 150 time points.  So the classifier is a movie.
%
beta = mdl.Beta;
nBasis = size(imageBasis,2);
img = zeros(cMosaic.rows,cMosaic.cols,tSamples);

colormap('default')
for ii=1:tSamples
    lst = (1:nBasis) + (ii-1)*nBasis;
    % lst = ii:tSamples:length(beta);
    tmp = imageBasis*beta(lst);
    img(:,:,ii) = reshape(tmp,cMosaic.rows,cMosaic.cols);
    imagesc(img(:,:,ii),[-.5 .5]); title(ii); colorbar; pause(0.2);
end
vcNewGraphWin; 
colormap('default'); ieMovie(img);

%%
