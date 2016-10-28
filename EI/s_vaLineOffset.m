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
ieInit;

%% Create display model
ppi    = 500;          % points per inch
imgFov = [.5 .5];      % image field of view
sensorFov = [.15 .15]; % field of view of retina
nFrames = 5000;        % Number of samples

vDist  = 0.3;          % viewing distance (meter)

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
tSamples = 100;  % 100 ms time sequence
prependZeros = 20;
weights = linspace(0, 1, tSamples/2 - prependZeros);
weights = [zeros(1, prependZeros) weights ...
           1-weights zeros(1, prependZeros)];

oiSeqAligned = oiSequence(OIs{3}, OIs{1}, 0.001, weights, ...
    'composition', 'blend');
oiSeqOffset = oiSequence(OIs{3}, OIs{2}, 0.001, weights, ...
    'composition', 'blend');

%%  Compute absorptions
nTrials = 1000;
dataAligned = [];
dataOffset = [];

cMosaic = coneMosaic;
cMosaic.setSizeToFOV(0.4 * imgFov);
cMosaic.integrationTime = 0.001;

for ii = 1 : nTrials
    cMosaic.emGenSequence(tSamples);
    emPos = cMosaic.emPositions;
    
    % compute for each time independently
    % will change to cMosaic.computeSeq
    curAligned = [];
    curOffset = [];
    for jj = 1 : tSamples
        cMosaic.emPositions = emPos(jj, :);
        cMosaic.compute(oiSeqAligned.frameAtIndex(jj), ...
            'currentFlag', false);
        curAligned = cat(3, curAligned, cMosaic.absorptions);
        cMosaic.compute(oiSeqOffset.frameAtIndex(jj), ...
            'currentFlag', false);
        curOffset = cat(3, curOffset, cMosaic.absorptions);
    end
    
    dataAligned = cat(1, dataAligned, curAligned(:)');
    dataOffset = cat(1, dataOffset, curOffset(:)');
    disp(['Trial: ' num2str(ii)]);
end

%% PCA and classification
% pca dimension reduction
nComponents = 50;
coefAligned = pca(dataAligned, 'NumComponents', nComponents);
pcaAligned = dataAligned * coefAligned;
coefOffset = pca(dataOffset, 'NumComponents', nComponents);
pcaOffset = dataOffset * coefOffset;

% svm classification
mdl = fitcsvm([pcaAligned; pcaOffset], ...
    [ones(nTrials, 1); -ones(nTrials, 1)], 'KernelFunction', 'linear');
crossMDL = crossval(mdl);
classLoss = kfoldLoss(crossMdl);