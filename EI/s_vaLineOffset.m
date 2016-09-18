%% s_vaLineOffset
%
%    Testing if people can see the difference between two cases:
%      1) A straight line
%      2) Two straight lines with 1 pixel apart
%    Vernier Acuity shows the positional acuity is around 6 sec of arc
%
% (HJ) Jan, 2014

%% Init Parameters
ieInit;

% Make a high resolution display for a small retinal field of view
ppi    = 500;          % points per inch
imgFov = [.5 .5];      % image field of view
sensorFov = [.15 .15]; % field of view of retina
nFrames = 5000;        % Number of samples

vDist  = 0.3;          % viewing distance (meter)

imgSz  = round(tand(imgFov)*vDist*(ppi/dpi2mperdot(1,'m')));   % number of pixels in image
imgFov = atand(max(imgSz)/ppi/(1/dpi2mperdot(1,'m'))/vDist);   % Actual fov

%% Create virtual display
display = displayCreate('LCD-Apple');
display = displaySet(display, 'dpi', ppi);

%% Create two scenes

% One with a constant line and one with an offset
scene = cell(2, 1);

params.display = display;
params.sceneSz = imgSz;
params.offset  = 0;
params.barWidth = 1;
params.bgColor = 0.5;
scene{1} = sceneCreate('vernier', 'display', params);

params.offset = 1;  % One pixel
scene{2} = sceneCreate('vernier', 'display', params);

% Uniform field ... we thinkg
params.barWidth = 0;
scene{3} = sceneCreate('vernier', 'display', params);

% set scene fov
for ii = 1 : 2
    scene{ii} = sceneSet(scene{ii}, 'h fov', imgFov);
    scene{ii} = sceneSet(scene{ii}, 'distance', vDist);
end

% Show radiance image (scene)
% vcAddObject(scene{1});
% vcAddObject(scene{2}); sceneWindow;

% Blank . Don't panic that the appearance is not gray.
% vcAddObject(scene{3}); sceneWindow;  

%% Create Human Lens

%  Create a typical human lens
% oi = oiCreate('wvf human',pupilMM,zCoefs,wave)
% Example:
%   To set the defocus, we can adjust one of the zCoefs (defocus).
%   We will rerun this with in focus and defocused case to check whether
%   this matters.
oi = oiCreate('wvf human');

% Compute optical image
OIs{1} = oiCompute(scene{1}, oi);
OIs{2} = oiCompute(scene{2}, oi);
OIs{3} = oiCompute(scene{3}, oi);

% Show irradiance (optical image) 
%vcAddObject(OIs{1}); oiWindow;
%vcAddObject(OIs{2}); oiWindow;
%vcAddObject(OIs{3}); oiWindow;

%%  Generate an eye movement sequence
cmosaic{1} = coneMosaic;
cmosaic{1}.emGenSequence(100);         % 
emPositions = cmosaic{1}.emPositions;

%% Noise free absorptions

% We create the absorptions with no noise and no eye movements
% We do this for a 1 ms exposure.  We will assemble the movies from these 1 ms
% noise free cases.
cmosaic{1} = coneMosaic;
cmosaic{2} = coneMosaic;
cmosaic{3} = coneMosaic;
LMS = cell(3,1);

for ii=1:3
    cmosaic{ii}.noiseFlag = false;  % There is no quantization in this case
    cmosaic{ii}.setSizeToFOV(imgFov);
    cmosaic{ii}.integrationTime = 0.001;   % ms
    
    % To account for eye movements, we need all of the L,M,S values at every
    % position.  We pull these out when we apply the emPath, below.
    LMS{ii} = cmosaic{ii}.computeSingleFrame(OIs{ii}, 'fullLMS', true);
end

%%  We might add a parameter to the window call to set a new figure
% cmosaic{1}.window;

%% Build up a movie that starts with uniform, puts on a stimulus, and ends with uniform

% The total movie is always 100 ms (frames)
% We put the stimulus on after 50 ms (frames)

% One trial

% We allow noise and eye movements
% Generate enough eye movements for all of the frames 

for ii= 1:3
    cmosaic{ii}.setSizeToFOV(sensorFov);
    
    % This creates a series of frames, shifted by the eye movements
    % Now add photon noise
    padRows = (size(LMS{1},1) - cmosaic{ii}.rows) / 2;
    padCols = (size(LMS{1},2) - cmosaic{ii}.cols) / 2;
    absorptions = cmosaic{ii}.applyEMPath(LMS{ii},'emPath',emPositions,...
        'padRows',padRows,'padCols',padCols);
    cmosaic{ii}.absorptions = cmosaic{ii}.photonNoise(absorptions);
end


%% Generate noise samples
%  Set exposure time to 1 ms
expTime = sensorGet(sensor, 'exp time');
emDuration = 0.001;
emPerExposure = expTime / emDuration;
sensor = sensorSet(sensor, 'exp time', emDuration);
sensor = sensorSetSizeToFOV(sensor, sensorFov(1), scene{1}, OIs{1});
sz = sensorGet(sensor, 'size');

% Generate eyemovement
p.nSamples = nFrames * 50;
sensor = eyemoveInit(sensor, p);

% Compute the cone absopritons
sensor = coneAbsorptions(sensor, OIs{1}, 2);

% Store the photon samples and add photons in one exposure time
pSamples1 = sensorGet(sensor, 'photons');
pSamples1 = sum(reshape(pSamples1, [sz nFrames 50]), 4);

% Compute cone absorptions for the second stimulus and store photon
% absorptions
sensor = coneAbsorptions(sensor, OIs{2}, 2);
pSamples2 = sensorGet(sensor, 'photons');
pSamples2 = sum(reshape(pSamples2, [sz nFrames 50]), 4);


%% Do it by SVM
% Classification
nFolds = 10;
labels = [ones(nFrames,1); -1*ones(nFrames,1)];
data = cat(1, RGB2XWFormat(pSamples1)', RGB2XWFormat(pSamples2)');

[acc, w] = svmClassifyAcc(data, labels, nFolds, 'linear');
fprintf('SVM acc:%f\n',acc(1));

% show weight image
vcNewGraphWin; imagesc(reshape(mean(w, 2), sz));

%%