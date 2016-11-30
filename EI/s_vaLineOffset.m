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

% Gaussian time series
tseries = ieScale(fspecial('gaussian',[1,150],30),0,1);

display = displayCreate('LCD-Apple');

clear sparams;  % Scene parameters in general
sparams.fov      = 0.35;  % Deg
sparams.distance = 0.3;   % Meters

% Basic vernier parameters for the oiSequence
clear vparams;
for ii = 3:-1:1
    vparams(ii) = vernierP;
    vparams(ii).display = display;
    vparams(ii).sceneSz =[50 50];  % This controls offset
end

% Uniform field
vparams(1).name = 'uniform'; vparams(1).bgColor = 0.5; vparams(1).barWidth = 0;

% Offset Line
vparams(2).name = 'offset';  vparams(2).bgColor = 0; vparams(2).offset = 1;

% Aligned lines
vparams(3).name = 'aligned'; vparams(3).bgColor = 0; vparams(3).offset = 0;

[offset, scenes] = oisCreate('vernier','add', tseries,'tparams',vparams([1 2]),'sparams',sparams);
offsetDeg = sceneGet(scenes{1},'degrees per sample')*vparams(2).offset;
fprintf('Offset in arc secs %.2f\n',offsetDeg*3600);

% offset.visualize;

aligned = oisCreate('vernier','add', tseries,'tparams',vparams([1 3]),'sparams',sparams);
% aligned.visualize;

%%  Compute absorptions
nTrials = 100;
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

% cMosaic.window;

%%  Reformat the time series for the PCA analysis

% imgListX matrix contains the temporal response for a pixel in a column
% The rows represent time samples by number of trials
% These are the temporal responses across all trials and time points.

rows = cMosaic.rows;
cols = cMosaic.cols;

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
tseries  = imgList * imageBasis;   %
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
nWeights = size(tseries,2);
data = zeros(2*nTrials,nWeights*tSamples);
for ii=1:(2*nTrials)
    start = (ii-1)*tSamples + 1;
    thisTrial = tseries(start:(start+tSamples - 1),:);
    data(ii,:) = thisTrial(:)';
end

mdl = fitcsvm(data,[ones(nTrials, 1); -ones(nTrials, 1)], 'KernelFunction', 'linear');

crossMDL = crossval(mdl);

func = @(y, yp, w, varargin) sum(abs(y(:, 1)-(yp(:, 1)>0)).*w);
classLoss = kfoldLoss(crossMDL, 'lossfun', func);
fprintf('Accuracy: %.2f%% \n', (1-classLoss) * 100);


%%
