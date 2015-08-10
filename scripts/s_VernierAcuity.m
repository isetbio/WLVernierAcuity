%% s_VernierAcuity
%    Testing if people can see the difference between two cases:
%      1) A straight line
%      2) Two straight lines with 1 pixel apart
%    Vernier Acuity shows the positional acuity is around 6 sec of arc
%
% (HJ) Jan, 2014


% Init Parameters
ieInit;

%% Simulation todo
%  1.  Run this with a monochrome eye for simplicity 
%  2.  Go back to 3 cones and change the stimulus SPD
%  3.  Change the eye movement parameters
%  4.  Change the WVF parameters for defocus
%  5.  Change the stimulus shape (e.g., use a harmonic)
%  6.  Illustrate the weights ...
%  7.  Change the SNR by changing the mean luminance level
%  8.  Change the contrast of the line
%  7.  Show the weights  as a function of the levels, parameters, so forth

tic

%% Set scene and sensor parameters
ppi = 1000;            % points per inch
imgFov = [.5 .5];      % image field of view
sensorFov = [.2 .2];   % field of view of sensor
nFrames = 5000;        % Number of samples

vDist  = 0.3;                                   % viewing distance (meter)
imgSz  = round(tand(imgFov)*vDist*39.37*ppi);   % number of pixels in image
imgFov = atand(max(imgSz)/ppi/39.37/vDist);     % Actual fov

%% Create virtual display
display = displayCreate('LCD-Apple');
display = displaySet(display, 'dpi', ppi);

%% Create vernier scenes

scene = cell(2, 1);

params.display = display;
params.sceneSz = imgSz;
params.offset  = 0;
params.barWidth = 1;
params.bgColor = 0.5;
scene{1} = sceneCreate('vernier', 'display', params);

params.offset = 1;
scene{2} = sceneCreate('vernier', 'display', params);

% set scene fov
for ii = 1 : 2
    scene{ii} = sceneSet(scene{ii}, 'h fov', imgFov);
    scene{ii} = sceneSet(scene{ii}, 'distance', vDist);
end

% Show radiance image (scene)
% vcAddAndSelectObject('scene', scene{1});
% vcAddAndSelectObject('scene', scene{2}); sceneWindow;

%% Create human optics

oi = oiCreate('wvf human');

% Compute optical image
OIs{1} = oiCompute(scene{1}, oi);
OIs{2} = oiCompute(scene{2}, oi);

% Show irradiance (optical image) 
%vcAddAndSelectObject('oi', OIs{1}); oiWindow;
%vcAddAndSelectObject('oi', OIs{2}); oiWindow;

%% Create Sensor
sensor = sensorCreate('human');

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
sensor1 = coneAbsorptions(sensor, OIs{1}, 2);

% Store the photon samples and add photons in one exposure time
pSamples1 = sensorGet(sensor1, 'photons');
pSamples1 = sum(reshape(pSamples1, [sz nFrames 50]), 4);

% Compute cone absorptions for the second stimulus and store photon
% absorptions
sensor2 = coneAbsorptions(sensor, OIs{2}, 2);
pSamples2 = sensorGet(sensor2, 'photons');
pSamples2 = sum(reshape(pSamples2, [sz nFrames 50]), 4);


%% SVM on the photon absorptions

% Classification parameters
nFolds = 10;
labels = [ones(nFrames,1); -1*ones(nFrames,1)];
data = cat(1, RGB2XWFormat(pSamples1)', RGB2XWFormat(pSamples2)');

% Calculate the svm.  This can take a while.
[acc, w] = svmClassifyAcc(data, labels, nFolds, 'linear');
fprintf('Cone photon absorptions (SVM acc):%f\n',acc(1));

% Show weight image
vcNewGraphWin; imagesc(reshape(mean(w, 2), sz));


%% SVM using the adapted cone current

% With cone noise added
params.addNoise = true;

[~,adapted1] = coneAdapt(sensor1,'rieke',params);
adapted1 = adapted1(:,:,1:50:end);

[~,adapted2] = coneAdapt(sensor2,'rieke',params);
adapted2 = adapted2(:,:,1:50:end);

% Classification folds
nFolds = 10;

% Spatial time by position
data = cat(1, RGB2XWFormat(adapted1)', RGB2XWFormat(adapted2)');

% One label for each time sample
labels = [ones(size(adapted1,3),1); -1*ones(size(adapted2,3),1)];

if ~isequal(length(labels), size(data,1))
    error('Something wrong in the label count');
end

% Perform the classification
[acc, w] = svmClassifyAcc(data, labels, nFolds, 'linear');

% Accuracy report
fprintf('Cone photocurrent (SVM acc):%f\n',acc(1));

vcNewGraphWin; imagesc(reshape(mean(w, 2), sz));

toc

%% Save critical parameters for reproducing this run and the output




%% END