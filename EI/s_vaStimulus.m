%% Template for making different OI sequences for vernier acuity
%
% The script explains the principles of oiSequence construction.
%
% BW, ISETBIO Team, 2016

%% Init Parameters
ieInit;

%% Create display model
ppi    = 500;          % points per inch
imgFov = [.5 .5];      % image field of view
sensorFov = [.15 .15]; % field of view of retina
% nFrames = 5000;        % Number of samples

vDist  = 0.3;          % viewing distance (meter)

% compute number of pixels in image
imgSz  = round(tand(imgFov)*vDist*(ppi/dpi2mperdot(1,'m')));

% compute actual fov of image, should be very close to desired fov
imgFov = atand(max(imgSz)/ppi/(1/dpi2mperdot(1,'m'))/vDist);

% Create virtual display
display = displayCreate('LCD-Apple', 'dpi', ppi);

%% Create scenes

% To create each oiSequence we need two types of scenes: one with a
% constant image and one with a time-varying presentation of the stimulus.
%
% We also want a pair of oiSequences, one for the offset and another with
% no offset.  We estimate the probability of discrimination by comparing
% these two sequences.
%

% Init scene cell array.  The background is constant and used for both of
% the oiSequences.  The two line scenes are added to the background.  One
% of the scenes has an offset and the other is a continuous line
%
% scene{1} - one aligned line
% scene{2} - two line segments with 1 pixel offset
% scene{3} - uniform background
scene = cell(3, 1);

% Init scene parameters
params.display = display;
params.sceneSz = imgSz;
params.barWidth = 1;

% Create scene{1}: constant line upon a background
params.offset  = 0;
params.bgColor = 0.0;
scene{1} = sceneCreate('vernier', 'display', params);
scene{1} = sceneSet(scene{1},'name','aligned');

% Create scene {2}: two line segments with 1 pixel offset on a background
params.offset = 1;  % units: pixel
params.bgColor = 0.0;
scene{2} = sceneCreate('vernier', 'display', params);
scene{2} = sceneSet(scene{2},'name','offset');

% Create scene {3}: background field only, no line
params.barWidth = 0;
params.bgColor = 0.5;
scene{3} = sceneCreate('vernier', 'display', params);
scene{3} = sceneSet(scene{3},'name','uniform');

% set scene fov
for ii = 1 : length(scene)
    scene{ii} = sceneSet(scene{ii}, 'h fov', imgFov);
    scene{ii} = sceneSet(scene{ii}, 'distance', vDist);
end

% Show radiance image (scene)
% vcAddObject(scene{1}); vcAddObject(scene{2}); sceneWindow; 

%% Compute human optical image sequences

% Create a typical human lens.  It will be possible to set the parameters
% for future experiments using
%
%   oi = oiCreate('wvf human',pupilMM,zCoefs,wave)
% 
% These are the default.
oi = oiCreate('wvf human');

% Compute optical images from the scene
OIs = cell(length(scene), 1);
for ii = 1 : length(OIs)
    OIs{ii} = oiCompute(scene{ii}, oi);
    OIs{ii} = oiCrop(OIs{ii},[8 8 52 52]);
end
% for ii=1:3, vcAddObject(OIs{ii}); end; oiWindow;

%% Build the aligned and offset oiSequences.

% We build the stimulus using a time series of weights. Then we make a
% linear ramp up for 30 ms, a ramp down for 30 ms, and we postpend zeros
% for 20 ms.
zTime = [50 150];
rTime = 15;
risingWeights = linspace(0, 1, rTime);   % 30 ms rising
weights = [zeros(1, zTime(1)) risingWeights ...
    1-risingWeights zeros(1, zTime(2))];
% vcNewGraphWin; plot(1:tSamples, weights,'o');
% xlabel('Time (ms)');

% Temporal samples.  Typically 1 ms, which is set by the parameter in the
% cone mosasic integration time.  That time is locked to the eye movements.
tSamples = length(weights); 
sampleTimes = 0.001*(1:tSamples);  % Time in sec

% The weights define some amount of the constant background and some amount
% of the line on the same constant background
oiSeqAligned = oiSequence(OIs{3}, OIs{1}, ...
    sampleTimes, weights, ...
    'composition', 'add');
oiSeqAligned.visualize('format','movie');

oiSeqOffset = oiSequence(OIs{3}, OIs{2}, ...
     sampleTimes, weights, ...
    'composition', 'add');
oiSeqOffset.visualize('format','movie');

%%



