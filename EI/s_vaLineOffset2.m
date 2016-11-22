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
vparams(1).barWidth = 0; vparams(1).bgColor = 0.1; vparams(1).name = 'uniform';

offset = oisCreate('vernier','add', weights,'hparams',vparams,'sparams',sparams);
offset.visualize;

vparams(2).name = 'aligned'; vparams(2).offset = 0;
aligned = oisCreate('vernier','add', weights,'hparams',vparams,'sparams',sparams);
aligned.visualize;

%%  Compute absorptions
nTrials = 3;
tSamples = aligned.length;

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
dataAligned  = cMosaic.compute(aligned,'emPaths',emPaths);
toc
% cMosaic.window;

%%  What PCA coefficient model should we use?

% This one simply strings out the time series points from all of the
% pixels, without any particular structure.  Every trial is then treated as
% the entries for the PCA analysis.  But the differences between the trials
% are mainly pixel-wise random.  So, we get very little savings.
% dataOffset  = reshape(dataOffset,nTrials,[]);
% dataAligned = reshape(dataAligned,nTrials,[]);

% This one asks whether we can replace the images by a smaller set of basis
% functions, and then we analyze the weights of these basis functions as a
% functon of time.
% Reshape dataAligned and dataOffset into matrices whos columns are images.
% Apply PCA, below.

% This is my effort.  Looks great!  The first 25 images look like the
% stimuli.

% Make a matrix with time x image (each row is an image).
vcNewGraphWin;
foo = zeros(nTrials*tSamples,37*37);
for tt = 1:nTrials
    lst = (1:tSamples) + (tSamples)*(tt-1);
    thisTrial = squeeze(dataAligned(tt,:,:,:));
    thisTrial = permute(thisTrial,[3 1 2]);
    thisTrial = reshape(thisTrial,tSamples,[]);
    foo(lst,:) = thisTrial;
end

% foo = reshape(dataAligned,nTrials,37*37,[]);
% foo = permute(foo,[1 3 2]);
% foo = reshape(foo,nTrials*110,[]);

% There are a couple of weird blanks.
% The ordering is not what I expect and doesn't correspond to the time
% series.
vcNewGraphWin;
for ii=size(foo,1)
    imagesc(reshape(foo(ii,:),37,37)); 
    drawnow; title(sprintf('%d',ii)); 
    pause(0.2); 
end

% Subtract the mean image from every row
foo = bsxfun(@minus,foo, mean(foo,1));

% Calculate the eigenvectors of the covariance matrix
[V,D] = eig(foo'*foo); 

% Show the eigenvalues.  The last values are the biggest.
D = diag(D);   
vcNewGraphWin; semilogy(D)

% Plot the largest eigenvectors as an image.
nComponents = 20;
for ii=1:nComponents
    colormap(gray(256));
    img = reshape(V(:,end-ii),37,37);
    imagesc(img);
    pause(0.2);
end

% Doesn't explain a lot of the variance.  But most of the variance is
% photon noise.  So, not sure what to think.
sum(D(end-nComponents:end)) / sum(D)

% The coefficient for the largest terms can be
%
%    coef = foo*V(:,largestNComponents);
%
% The matrix coef will be time x nComponents
% The reconstruction of the original image will be
%
%   reconstructed  = V(:largestNComponents)*coef'; 
%                  = coef*V(:,largestNComponents)';  % time x nPixels
%
%   reshape(reconstructed)
%
% Where reconstructed will be time x


%%
tic
dataOffset = cMosaic.compute(offset,'emPaths',emPaths);
toc
% cMosaic.window;

foo = reshape(dataOffset,100,37*37,[]);
foo = permute(foo,[1 3 2]);
foo = reshape(foo,100*110,[]);
nComponents = 10;
[coefOffset,~,~,~,explainedOffset] = pca(foo, 'NumComponents', nComponents);

foo = bsxfun(@minus,foo, mean(foo,1));
[V,D] = eig(foo'*foo); 
D = diag(D);   % The last values are the biggest.

vcNewGraphWin; semilogy(D)

for ii=1:25
    colormap(gray(256));
    img = reshape(V(:,end-ii),37,37);
    imagesc(img);
    pause(1);
end

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
