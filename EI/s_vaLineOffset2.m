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

%% Init Parameters

stimWeights = ieScale(fspecial('gaussian',[1,50],15),0,1);
weights = [zeros(1, 30), stimWeights, zeros(1, 30)];

clear sparams;
sparams.fov = 0.5;

clear vparams; 
vparams(2) = vernierP;
vparams(2).name = 'offset'; vparams(2).bgColor = 0; vparams(2).offset = 5;
vparams(1) = vparams(2);
vparams(1).barWidth = 0; vparams(1).bgColor = 0.5; vparams(1).name = 'uniform';

offset = oisCreate('vernier','add', weights,'hparams',vparams,'sparams',sparams);
offset.visualize;

vparams(2).name = 'aligned'; vparams(2).offset = 0;
aligned = oisCreate('vernier','add', weights,'hparams',vparams,'sparams',sparams);
aligned.visualize;

%%  Compute absorptions
nTrials = 100;

cMosaic = coneMosaic;
cMosaic.setSizeToFOV(0.4 * oiGet(aligned.oiFixed,'fov'));
cMosaic.integrationTime = 0.001;

emPaths = zeros(nTrials, tSamples, 2);
for ii = 1 : nTrials
    cMosaic.emGenSequence(tSamples);
    emPaths(ii, :, :) = cMosaic.emPositions;
end

%%

tic
[dataOffset]  = cMosaic.compute(aligned,'emPaths',emPaths);
toc
% cMosaic.window;

tic
[dataAligned] = cMosaic.compute(offset,'emPaths',emPaths);
toc
% cMosaic.window;


%%  What PCA coefficient model should we use?

% This one simply strings out the time series points from all of the
% pixels, without any particular structure.  Every trial is then treated as
% the entries for the PCA analysis.  But the differences between the trials
% are mainly pixel-wise random.  So, we get very little savings.
dataOffset  = reshape(dataOffset,nTrials,[]);
dataAligned = reshape(dataAligned,nTrials,[]);

% This one asks whether we can replace the images by a smaller set of basis
% functions, and then we analyze the weights of these basis functions as a
% functon of time.
% Reshape dataAligned and dataOffset into matrices whos columns are images.
% Apply PCA, below.

%%

fprintf('Dimension reduction with PCA...');
nComponents = 30;
[coefAligned,~,~,~,explainedAligned] = pca(dataAligned, 'NumComponents', nComponents);
pcaAligned  = dataAligned * coefAligned;
sum(explainedAligned(1:30))

coefOffset  = pca(dataOffset, 'NumComponents', nComponents);
pcaOffset   = dataOffset * coefOffset;
fprintf('Done\n');

% If you want to visualize some ...
%
% plot3(pcaOffset(:,1),pcaOffset(:,2),pcaOffset(:,3),'ro', ...
%     pcaAligned(:,1),pcaAligned(:,2),pcaAligned(:,3),'bo');
% grid on

%% svm classification
fprintf('SVM Classification ');
mdl = fitcsvm([pcaAligned; pcaOffset], ...
    [ones(nTrials, 1); -ones(nTrials, 1)], 'KernelFunction', 'linear');
crossMDL = crossval(mdl);
classLoss = kfoldLoss(crossMDL);
fprintf('Accuracy: %.2f%% \n', (1-classLoss) * 100);

%%
