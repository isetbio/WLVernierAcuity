%% s_vaAbsorptions
%
%    Testing if people can see the difference between two cases:
%      1) A straight line
%      2) Two straight lines barOffset pixels apart
%
%  Vernier acuity in human shows the positional acuity is around 6 sec of
%  arc. Here, we are analyzing stimulus, optics, and eye movements and
%  basing the calculation on absorptions.
%
%  In a separate script, we analyze the discriminability based on
%  photocurrent rather than absorptions.
%
% In this case we try
%    Standard retinal parameters
%    White line on a gray monitor background
%    Sweep out viewing distance
%
% HJ/BW, ISETBIO TEAM, 2016

%% Moved parameter setting to s_EIScratch.m
%
% We now run this script just to execute.
%

%% If already computed, use the imageBasis.  If not, make an image basis.
tmp = [];
if exist('imageBasisAbsorptions.mat','file'), tmp = load('imageBasisAbsorptions'); end
if isfield(tmp,'basisParameters')
    basisParameters = tmp.basisParameters;
    if      isequal(basisParameters.timesd, params.timesd) && ...
            isequal(basisParameters.tStep,params.tStep) && ...
            isequal(basisParameters.fov,params.fov) && ...
            isequal(basisParameters.vernier.bgColor, params.vernier.bgColor) && ...
            isequal(basisParameters.vernier.barWidth,params.vernier.barWidth) && ...
            isequal(basisParameters.vernier.barLength,params.vernier.barLength) && ...
            isequal(basisParameters.vernier.gap,params.vernier.gap) && ...
            isequal(basisParameters.vernier.barColor,params.vernier.barColor) && ...
            isequal(basisParameters.vernier.sceneSz,params.vernier.sceneSz)

        disp('Loading image basis because parameters match')
        load('imageBasisAbsorptions','imageBasis');
    else
        disp('Creating new image basis - parameters do not match')
        imageBasis = vaPCA(params);
    end
else
    disp('Creating new image basis - can not find parameters in file')
    imageBasis = vaPCA(params);
end
imageBasis = imageBasis(:,1:params.nBasis);

% % % Have a look if you like
% vcNewGraphWin; colormap(gray(256))
% mx = max(imageBasis(:)); mn = min(imageBasis(:));
% for ii=1:params.nBasis
%     imagesc(reshape(imageBasis(:,ii),cMosaic.rows,cMosaic.cols),[mn mx]);
%     title(sprintf('Basis %d',ii));
%     pause(0.5);
% end


%% Create the aligned and offset vernier stimuli

% This could loop here on the barOffset
X = zeros(1,numel(barOffset));
P = zeros(1,numel(barOffset));

%%
for bb = 1:numel(barOffset)
    
    params.vernier.offset = barOffset(bb);
    [aligned, offset, scenes,tseries] = vaStimuli(params);
    % offset.visualize;
    % aligned.visualize;
    % ieAddObject(offset.oiModulated); oiWindow;
    % ieAddObject(scenes{2}); sceneWindow;
    % vcNewGraphWin; plot(params.tsamples,tseries)
    
    % Offset lines
    % offsetDeg = sceneGet(scenes{1},'degrees per sample')*vparams(2).offset;
    % fprintf('Offset in arc secs %.2f\n',offsetDeg*3600);
    
    %%  Compute absorptions for multiple trials
    
    tSamples = aligned.length;
    
    cMosaic = coneMosaic;
    
    % Set the mosaic size to 15 minutes (.25 deg) because that is the spatial
    % pooling size found by Westheimer and McKee
    cMosaic.setSizeToFOV(params.fov);
    
    % Not sure why these have to match, but there is a bug if they don't.
    cMosaic.integrationTime = aligned.timeStep;
    
    %% For aligned or offset
    
    cMosaic.noiseFlag    = 'random';
    
    emPaths  = cMosaic.emGenSequence(tSamples,'nTrials',nTrials,'em',params.em);
    
    disp('Aligned')
    tic
    alignedA = cMosaic.compute(aligned,'currentFlag',false,'emPaths',emPaths);
    toc
    
    disp('Offset')
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
    
    % Could shrink nBasis here ... or not
    % % imageBasis = imageBasis(:,1:nBasis);
    
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
    
    % func = @(y, yp, w, varargin) sum(abs(y(:, 1)-(yp(:, 1)>0)).*w);
    % classLoss = kfoldLoss(crossMDL, 'lossfun', func);
    
    % Select some of the data (80%) as the training set.
    train_index = zeros(nTrials, 1);
    train_index(randperm(nTrials, round(0.8*nTrials))) = 1;
    train_index = train_index > 0;
    
    % The aligned and offset trials are still matched
    train_index = [train_index; train_index];
    
    % Fit the SVM model.
    mdl = fitcsvm(data(train_index, :), label(train_index), ...
        'KernelFunction', 'linear');
    
    % predict the data not in the training set.
    yp = predict(mdl, data(~train_index, :));
    classLoss = sum(label(~train_index) ~= yp) / length(yp);
    X(bb) = barOffset(bb); P(bb) = (1-classLoss) * 100;
    fprintf('Accuracy for held out data: (%d, %.2f%%) \n', barOffset(bb), (1-classLoss) * 100);
end

%% Dump out for saving and plotting

params
s = sprintf('X = ['); s = [s, sprintf('%d ',X)]; s = [s , sprintf(']')];
disp(s)
s = sprintf('P = ['); s = [s, sprintf('%.2f ',P)]; s = [s , sprintf(']')];
disp(s)

vcNewGraphWin;
plot(6*X,P,'o-');
grid on
xlabel('Offset (arc sec)'); ylabel('Percent correct');
set(gca,'ylim',[45 100])
title(sprintf('Bar width %.1f (sec), duration %.3f (sec) (A)',6*params.vernier.barWidth, params.timesd));

%% Run cross validation

% I don't understand this yet, so I am skipping (BW)

% crossMDL = crossval(mdl);    % I don't understand this syntax
% func = @(y, yp, w, varargin) sum(abs(y(:, 1)-(yp(:, 1)>0)).*w);
% classLoss = kfoldLoss(crossMDL, 'lossfun', func);

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
%  When we get mdl.Beta, it has size of tSamples*nBasis.   So, if
%  there are 150 tSamples and 10 spatial image components, then Beta is
%  1500.
%
%  I think that on a single trial we have tSamples*nBasis numbers.  If we
%  want to express these  numbers as a movie, we multiply the weights times
%  the spatial image basis at each of the time points. So the classifier is
%  a movie.
%
%
%
% % These are the beta weights we multiply
% beta = mdl.Beta;
% % nBasis = size(imageBasis,2);
% img = zeros(cMosaic.rows,cMosaic.cols,tSamples);
%
% colormap('default')
% for ii=1:tSamples
%     lst = (1:nBasis) + (ii-1)*nBasis;
%     % lst = ii:tSamples:length(beta);
%     tmp = imageBasis*beta(lst);
%     img(:,:,ii) = reshape(tmp,cMosaic.rows,cMosaic.cols);
%     imagesc(img(:,:,ii),[-.5 .5]); title(ii); colorbar; pause(.2);
% end
%
% vcNewGraphWin;
% colormap('default'); ieMovie(img);

%%