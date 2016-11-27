%% s_vaLineOffset2 - Will become LineOffset again.  Just simplifying.
%
%    Testing if people can see the difference between two cases:
%      1) A straight line
%      2) Two straight lines with 1 pixel apart
%    Vernier Acuity shows the positional acuity is around 6 sec of arc
%
% In this case we try
%
%    Standard retinal parameters
%    White line on a gray monitor background
%    Sweep out viewing distance
%
% We are running in the ConeMosaicOSintegrationTime branch, which will
% probably become the master branch before too long
%
% HJ/BW, ISETBIO TEAM, 2016

%%
ieInit

%% Init Parameters

% Think about the time series here.
stimWeights = ieScale(fspecial('gaussian',[1,50],15),0,1);
weights = [zeros(1, 30), stimWeights, zeros(1, 30)];

clear sparams;
sparams.fov      = 0.35;

clear vparams; 
vparams(2) = vernierP;
vparams(2).name = 'offset'; vparams(2).bgColor = 0; vparams(2).offset = 1;
vparams(1) = vparams(2);
vparams(1).barWidth = 0; vparams(1).bgColor = 0.5; vparams(1).name = 'uniform';

offset = oisCreate('vernier','add', weights,...
    'tparams',vparams,'sparams',sparams);
% offset.visualize;

vparams(2).name = 'aligned'; vparams(2).offset = 0;
aligned = oisCreate('vernier','add', weights,...
    'tparams',vparams,'sparams',sparams);
% aligned.visualize;

%%  Compute absorptions
nTrials = 200;
tSamples = aligned.length;

cMosaic = coneMosaic;
cMosaic.setSizeToFOV(0.4 * oiGet(aligned.oiFixed,'fov'));
cMosaic.integrationTime = 0.001;

%% For aligned or offset

tic

% Eye movements and computations for each type of stimulus
emPaths = zeros(nTrials, tSamples, 2);
for ii = 1 : nTrials
    cMosaic.emGenSequence(tSamples);
    emPaths(ii, :, :) = cMosaic.emPositions;
end

cMosaic.os.noiseFlag = true;
alignedA = cMosaic.compute(aligned, ...
    'emPaths',emPaths, ...
    'currentFlag',false);

% Separate eye movements
emPaths = zeros(nTrials, tSamples, 2);
for ii = 1 : nTrials
    cMosaic.emGenSequence(tSamples);
    emPaths(ii, :, :) = cMosaic.emPositions;
end

cMosaic.os.noiseFlag = true;
offsetA = cMosaic.compute(offset, ...
    'emPaths',emPaths, ...
    'currentFlag',false);
toc
rows = cMosaic.rows;
cols = cMosaic.cols;

% cMosaic.window;

%%  What PCA coefficient model should we use?

% We replace the images by a smaller set of basis functions. Then we
% analyze the weights of these basis functions as a function of time.
%
% BW should comment more.
% outData is trial x X x Y x T
% We want (trial*Time) x (space)
%
% Make a matrix with time x image (each row is an image).
imgListAligned = zeros(nTrials*tSamples,rows*cols);
for tt = 1:nTrials
    lst = (1:tSamples) + (tSamples)*(tt-1);
    thisTrial = squeeze(alignedA(tt,:,:,:));
    thisTrial = permute(thisTrial,[3 1 2]);  
    thisTrial = reshape(thisTrial,tSamples,[]);
    imgListAligned(lst,:) = thisTrial;
end

imgListOffset = zeros(nTrials*tSamples,rows*cols);
for tt = 1:nTrials
    lst = (1:tSamples) + (tSamples)*(tt-1);
    thisTrial = squeeze(offsetA(tt,:,:,:));
    thisTrial = permute(thisTrial,[3 1 2]);  
    thisTrial = reshape(thisTrial,tSamples,[]);
    imgListOffset(lst,:) = thisTrial;
end

% Visualize
%
% imgList = abs(imgListOffset);  % Absorptions or current should be positive
% mx = max(imgListOffset(:));
% vcNewGraphWin; colormap(gray(mx));
% for ii=1:size(imgListOffset,1)
%     image(reshape(imgListOffset(ii,:),rows*cols)); 
%     drawnow; title(sprintf('%d',ii)); 
%     pause(0.05); 
% end

%% Home grown not-centered PCA

% We make bases for each type of stimulus and then concatenate them.
nComponents = 3;
[~,~,V] = svd(imgListAligned,'econ');
imageBasis = V(:,1:nComponents);

[~,~,V] = svd(imgListOffset,'econ');
imageBasis = cat(2,imageBasis,V(:,1:nComponents));

% vcNewGraphWin; colormap(gray(256))
% for ii=1:(2*nComponents)
%     imagesc(reshape(imageBasis(:,ii),rows*cols));
%     pause(0.5);
% end

imgList = cat(1,imgListAligned,imgListOffset);
weights  = imgList * imageBasis;   %
% plot(weights);

%% Let's reconstruct the approximate absorption sequence
% 
% recon = imageBasis*weights';
% vcNewGraphWin; colormap(gray(round(max(recon(:)))));
% for ii=1:size(recon,2)
%     image(reshape(recon(:,ii),rows*cols));
%     title(sprintf('%d',ii));
%     pause(0.1);
% end


%% svm classification

fprintf('SVM Classification ');

% Put the weights from each trial into the rows of a matrix
% Each row is another trial
nWeights = size(weights,2);
data = zeros(2*nTrials,nWeights*tSamples);
for ii=1:(2*nTrials)
    start = (ii-1)*tSamples + 1;
    thisTrial = weights(start:(start+tSamples - 1),:);
    data(ii,:) = thisTrial(:)';
end

mdl = fitcsvm(data,[ones(nTrials, 1); -ones(nTrials, 1)], 'KernelFunction', 'linear');

crossMDL = crossval(mdl);

func = @(y, yp, w, varargin) sum(abs(y(:, 1)-(yp(:, 1)>0)).*w);
classLoss = kfoldLoss(crossMDL, 'lossfun', func);
fprintf('Accuracy: %.2f%% \n', (1-classLoss) * 100);


%%
