%% t_VernierClassifier
%
% This tutorial script uses biological and computational methods to analyze
% vernier acuity (super acuity) in human vision.
%
% This script will explain how we use SVM classifier to partially explain
% the hyper-acuity story.
%
%  HJ/BW, ISETBIO TEAM, 2015

%  Initialize a new session
ieInit;

%% Create the display
% In this example we impose a linear gamma table, though
% in general it could be the default or anything.
dpi = 500; d = displayCreate('LCD-Apple','dpi',dpi);

viewDist = 2; % viewing distance in meters
d = displaySet(d, 'viewing distance', viewDist);
d = displaySet(d, 'gamma', 'linear');

%% Create Vernier Scene (full display radiance representation)
[~, p] = imageVernier();   % Mainly to get the parameters
p.pattern = 0.2*ones(1,513); p.pattern(257) = 1;
p.sceneSz = [513 513];

% Aligned
p.offset = 0;
imgA = imageVernier(p);

% Misaligned
p.offset = 2;
imgM = imageVernier(p);
        
% Create a scene with the image using the display parameters
% The scene spectral radiance is created using the RGB image and the
% properties of the display.
sceneA = sceneFromFile(imgA, 'rgb', [], d); % aligned
sceneM = sceneFromFile(imgM, 'rgb', [], d); % mis-aligned

fov = size(imgA,2)/displayGet(d,'dots per deg');
sceneA = sceneSet(sceneA,'fov',fov);
sceneM = sceneSet(sceneM,'fov',fov);

%% Compute Irradiance with Optics Wavefront
oi = oiCreate('wvf human');
oiA = oiCompute(sceneA, oi);
oiM = oiCompute(sceneM, oi);

%% Compute Photon Absorptions of Human Cone Receptors
% We compute a lot of 50 ms cone absorptions here
% And, it's equivalent to a series of experiment with 50 ms stimulus
% in each trial
n_samples = 3000;

% Compute the human cone absorption samples no eye movement.
sensor = sensorCreate('human');
fov = 0.4 * sceneGet(sceneA, 'fov');
sensor = sensorSetSizeToFOV(sensor, fov, sceneA, oiA);
sensor = sensorSet(sensor, 'positions', zeros(n_samples, 2));

sensor = coneAbsorptions(sensor, oiA);
photons_A = sensorGet(sensor, 'photons');

sensor = coneAbsorptions(sensor, oiM);
photons_M = sensorGet(sensor, 'photons');

%% Classify based on no eye-movement cone absorptions
% set up labels for two groups, 1 for aligned, -1 for mis-aligned
labels = [ones(n_samples, 1); -1*ones(n_samples, 1)];

% prepare photon absorption data
% set cone absorptions in one integration time in rows of data matrix
data = [RGB2XWFormat(photons_A)'; RGB2XWFormat(photons_M)'];

% find a linear boundray that best split the two groups with svm
% here, we do a 10 fold cross-validation
nFolds = 10;
[acc1, w1] = svmClassifyAcc(data, labels, nFolds, 'linear');

% let's have a look at the weights
vcNewGraphWin; imagesc(reshape(mean(w1, 2), sensorGet(sensor, 'size')));

% Without eye movement, only the position of the top bar is used
% However, because of the eye movement, we cannot perfect fix our eyes on
% the desired point. The actual position of the eye when stimulus first
% shown is randomly distributed around the fixation point

%% Randomness of initial eye position
%  now suppose the eye is at random position at the begin and hold that
%  position on for the entire 50 ms.
sensor = sensorSet(sensor, 'positions', round(2*randn(n_samples, 2)));

sensor = coneAbsorptions(sensor, oiA);
photons_A = sensorGet(sensor, 'photons');

sensor = coneAbsorptions(sensor, oiM);
photons_M = sensorGet(sensor, 'photons');

data = [RGB2XWFormat(photons_A)'; RGB2XWFormat(photons_M)'];
[acc2, w2] = svmClassifyAcc(data, labels, nFolds, 'linear');

% let's have a look at the weights again
vcNewGraphWin; imagesc(reshape(mean(w2, 2), sensorGet(sensor, 'size')));

%% Using real fixational eye movement
%  init eye movement path
exp_time = sensorGet(sensor, 'exp time');
sample_time = 0.001; % we always use 1ms to generate eye movement
p.nSamples = n_samples * exp_time / sample_time;
sensor = eyemoveInit(sensor, p);
sensor = sensorSet(sensor, 'exp time', sample_time);

%  compute cone absorption
sensor = coneAbsorptions(sensor, oiA);
photons_A = sensorGet(sensor, 'photons');

sensor = coneAbsorptions(sensor, oiM);
photons_M = sensorGet(sensor, 'photons');

%  now we have data integrated by every 1 ms, integrate to 50 ms
photons_A = sum(reshape(photons_A, [sz n_samples exp_time/sample_time]), 4);
photons_M = sum(reshape(photons_M, [sz n_samples exp_time/sample_time]), 4);

%  do classification again
data = [RGB2XWFormat(photons_A)'; RGB2XWFormat(photons_M)'];
[acc3, w3] = svmClassifyAcc(data, labels, nFolds, 'linear');

% now the accuracy is lower
vcNewGraphWin; imagesc(reshape(mean(w2, 2), sensorGet(sensor, 'size')));