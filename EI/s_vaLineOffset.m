%% s_vaLineOffset
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
% ieInit;
if notDefined('computeStage')
    computeStage = 'absorption';  % 'absorption', 'current'
end
if strcmp(computeStage, 'current')
    error('Computing current for each time (implemented below) is wrong');
end    

%% Create display model
ppi    = 500;          % points per inch
imgFov = [.5 .5];      % image field of view
sensorFov = [.15 .15]; % field of view of retina
nFrames = 5000;        % Number of samples

if ~exist('vDist', 'var')
vDist  = 0.3;          % viewing distance (meter)
end

% compute number of pixels in image
imgSz  = round(tand(imgFov)*vDist*(ppi/dpi2mperdot(1,'m')));

% compute actual fov of image, should be very close to desired fov
imgFov = atand(max(imgSz)/ppi/(1/dpi2mperdot(1,'m'))/vDist);

% Create virtual display
display = displayCreate('LCD-Apple', 'dpi', ppi);

%% Create scenes
% In general, we want two types of scenes: one with a constant line 
% and one with an offset. In order to handle temporal variation, we create
% a background scene, which is a uniform field. We use the the background
% to blend with the two scenes and get the temporal varying OI sequence.

% Init scene array
% scene{1} - constant line
% scene{2} - two line segments with 1 pixel offset
% scene{3} - uniform background
scene = cell(3, 1);

% Init scene parameters
params.display = display;
params.sceneSz = imgSz;
params.offset  = 0;
params.barWidth = 1;
params.bgColor = 0.5;

% Create scene{1}: constant line
scene{1} = sceneCreate('vernier', 'display', params);

% Create scene {2}: two line segments with 1 pixel offset
params.offset = 1;  % units: pixel
scene{2} = sceneCreate('vernier', 'display', params);

% Create scene {3}: uniform field
params.barWidth = 0;
scene{3} = sceneCreate('vernier', 'display', params);

% set scene fov
for ii = 1 : length(scene)
    scene{ii} = sceneSet(scene{ii}, 'h fov', imgFov);
    scene{ii} = sceneSet(scene{ii}, 'distance', vDist);
end

% Show radiance image (scene)
% vcAddObject(scene{1}); vcAddObject(scene{2}); sceneWindow; 

%% Compute human optical image
% Create a typical human lens
%   oi = oiCreate('wvf human',pupilMM,zCoefs,wave)
oi = oiCreate('wvf human');

% Compute optical image
OIs = cell(length(scene), 1);
for ii = 1 : length(OIs)
    OIs{ii} = oiCompute(scene{ii}, oi);
end

% Build oiSequence
tSamples = 50;  % 100 ms time sequence
prependZeros = 10;
weights = linspace(0, 1, tSamples/2 - prependZeros);
weights = [zeros(1, prependZeros) weights ...
           1-weights zeros(1, prependZeros)];

oiSeqAligned = oiSequence(OIs{3}, OIs{1}, 0.001*(1:tSamples), weights, ...
    'composition', 'blend');
oiSeqOffset = oiSequence(OIs{3}, OIs{2}, 0.001*(1:tSamples), weights, ...
    'composition', 'blend');

%%  Compute absorptions
nTrials = 200;
dataAligned = [];
dataOffset = [];

cMosaic = coneMosaic;
cMosaic.setSizeToFOV(0.4 * imgFov);
cMosaic.integrationTime = 0.001;

emPos = zeros(tSamples, 2, nTrials);
for ii = 1 : nTrials
    cMosaic.emGenSequence(tSamples);
    emPos(:, :, ii) = cMosaic.emPositions;
end

fprintf('Computing cone response: ');
msg = '';
for ii = 1 : tSamples
    % compute for all trials at this time point
    cMosaic.emPositions = squeeze(emPos(ii, :, :))';
    if strcmp(computeStage, 'absorption')
        % compute photon absorptions at time ii
        cMosaic.compute(oiSeqAligned.frameAtIndex(ii), ...
            'currentFlag', false);
        dataAligned = cat(2, dataAligned, ...
            reshape(cMosaic.absorptions, [], nTrials)');
        cMosaic.compute(oiSeqOffset.frameAtIndex(ii), ...
            'currentFlag', false);
        dataOffset = cat(2, dataOffset, ...
            reshape(cMosaic.absorptions, [], nTrials)');
    elseif strcmp(computeStage, 'current')
        % compute photo current at time ii
        cMosaic.compute(oiSeqAligned.frameAtIndex(ii));
        dataAligned = cat(2, dataAligned, ...
            reshape(cMosaic.current, [], nTrials)');
        cMosaic.compute(oiSeqOffset.frameAtIndex(ii));
        dataOffset = cat(2, dataOffset, ...
            reshape(cMosaic.current, [], nTrials)');
    end     
    
    fprintf(repmat('\b', [1 length(msg)]));
    msg = sprintf('%d / %d', ii, tSamples);
    fprintf(msg);
end
fprintf('\n');

%% PCA and classification
% pca dimension reduction

fprintf('Dimension reduction with PCA...');
nComponents = 20;
coefAligned = pca(dataAligned, 'NumComponents', nComponents);
pcaAligned = dataAligned * coefAligned;
coefOffset = pca(dataOffset, 'NumComponents', nComponents);
pcaOffset = dataOffset * coefOffset;
fprintf('Done\n');

% svm classification
fprintf('SVM Classification ');
mdl = fitcsvm([pcaAligned; pcaOffset], ...
    [ones(nTrials, 1); -ones(nTrials, 1)], 'KernelFunction', 'linear');
crossMDL = crossval(mdl);

func = @(y, yp, w, varargin) sum(abs(y(:, 1)-(yp(:, 1)>0)).*w);
classLoss = kfoldLoss(crossMDL, 'lossfun', func);
fprintf('Accuracy: %.2f%% \n', (1-classLoss) * 100);


