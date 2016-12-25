%% s_vaAbsorptions
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
nTrials = 100;
tStep   = 5;

% Set basic parameters for the vernier stimulus
clear p; 
p.barOffset = 0;     % Pixels on the display
p.barWidth  = 1;     % Pixels on the display
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

% Not sure why these have to match, but there is a bug if they don't.
cMosaic.integrationTime = aligned.timeStep;  

%% For aligned or offset

cMosaic.noiseFlag    = 'random';

emPaths  = cMosaic.emGenSequence(tSamples,'nTrials',nTrials);

tic
alignedA = cMosaic.compute(aligned,'currentFlag',false,'emPaths',emPaths);
toc

tic
offsetA = cMosaic.compute(offset,'currentFlag',false,'emPaths',emPaths);
toc

% cMosaic.window;

%%  Reformat the time series for the PCA analysis

% imgListX matrix contains the temporal response for a pixel in a column
% The rows represent time samples by number of trials
% These are the temporal responses across all trials and time points.

imgListAligned = trial2Matrix(alignedA,cMosaic);
imgListOffset  = trial2Matrix(offsetA,cMosaic);

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

% % We make bases for each type of stimulus and then concatenate them.
% nComponents = 5;
% [~,~,V] = svd(imgListAligned,'econ');
% imageBasis = V(:,1:nComponents);
% 
% [~,~,V] = svd(imgListOffset,'econ');
% imageBasis = cat(2,imageBasis,V(:,1:nComponents));
% 
% % Have a look if you like
% % vcNewGraphWin; colormap(gray(256))
% % for ii=1:(2*nComponents)
% %     imagesc(reshape(imageBasis(:,ii),rows*cols));
% %     pause(0.5);
% % end
% 

% produced by s_vaAbsorptionsPCA.  Needs to match the stimulus here, I
% think.
load('imageBasisAbsorptions','imageBasis');

imgList = cat(1,imgListAligned,imgListOffset);
% imgList = cat(1,imgListAligned,imgListAligned);

% Time series of weights
weightSeries  = imgList * imageBasis;  

%% Let's see if we can reduce the dimensionality of the time series 
%
% The reasons is that the photocurrent time series, which smooths the
% signal over time, performs much better with the SVM.  So, I think that
% smoothing the time series in the absorptions would allow the SVM to find
% a good solution here, as well.
%
% Now the weight series for each image basis comprises 150 numbers for each
% of 600 trials. We frame this as 150 x 600 and reduce it to [150 x
% nTimeBasis]*wgts

% These are the time series for each of the trials
% foo = reshape(weightSeries(:,1),tSamples,2*nTrials);
% vcNewGraphWin; plot(foo);
% [U,S,T] = svd(foo,'econ');
% vcNewGraphWin;
% plot(U(:,1));
% 
% wgt = U(:,1:3)'*foo;
% vcNewGraphWin; 
% plot3(wgt(1,1:300),wgt(2,1:300),wgt(3,1:300),'ro')
% hold on; 
% plot3(wgt(1,301:600),wgt(2,301:600),wgt(3,301:600),'go')
% hold off

%% Start classification training

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
label = [ones(nTrials, 1); -ones(nTrials, 1)];

%
% func = @(y, yp, w, varargin) sum(abs(y(:, 1)-(yp(:, 1)>0)).*w);
% classLoss = kfoldLoss(crossMDL, 'lossfun', func);
train_index = zeros(nTrials, 1);
train_index(randperm(nTrials, 0.8*nTrials)) = 1;
train_index = train_index > 0;
train_index = [train_index; train_index];

mdl = fitcsvm(data(train_index, :), label(train_index), ...
    'KernelFunction', 'linear');

% predict on training set - I think we should generate a completely fresh
% test here. (BW)
yp = predict(mdl, data(~train_index, :));
classLoss = sum(label(~train_index) ~= yp) / length(yp);
fprintf('Accuracy: %.2f%% \n', (1-classLoss) * 100);

% cross validation
crossMDL = crossval(mdl);
func = @(y, yp, w, varargin) sum(abs(y(:, 1)-(yp(:, 1)>0)).*w);
classLoss = kfoldLoss(crossMDL, 'lossfun', func);

%%
