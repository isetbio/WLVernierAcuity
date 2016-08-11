%% s_VernierAcuity
%    Testing if people can see the difference between two cases:
%      1) A straight line
%      2) Two straight lines with 1 pixel apart
%    Vernier Acuity shows the positional acuity is around 6 sec of arc
%
% (HJ) Jan, 2014

%% Init Parameters
ieInit;

ppi    = 500;          % points per inch
imgFov = [.5 .5];      % image field of view
sensorFov = [.15 .15]; % field of view of sensor
nFrames = 5000;        % Number of samples

vDist  = 0.3;                                   % viewing distance (meter)

imgSz  = round(tand(imgFov)*vDist*(ppi/dpi2mperdot(1,'m')));   % number of pixels in image
imgFov = atand(max(imgSz)/ppi/(1/dpi2mperdot(1,'m'))/vDist);    % Actual fov

%% Create virtual display
display = displayCreate('LCD-Apple');
display = displaySet(display, 'dpi', ppi);

%% Create Scene
%  create vernier scene
scene = cell(2, 1);

params.display = display;
params.sceneSz = imgSz;
params.offset  = 0;
params.barWidth = 1;
params.bgColor = 0.5;
scene{1} = sceneCreate('vernier', 'display', params);

params.offset = 1;  % One pixel
scene{2} = sceneCreate('vernier', 'display', params);

% set scene fov
for ii = 1 : 2
    scene{ii} = sceneSet(scene{ii}, 'h fov', imgFov);
    scene{ii} = sceneSet(scene{ii}, 'distance', vDist);
end

% Show radiance image (scene)
% vcAddObject(scene{1});
% vcAddObject(scene{2}); sceneWindow;

%% Create Human Lens
%  Create a typical human lens
oi = oiCreate('wvf human');

% Compute optical image
OIs{1} = oiCompute(scene{1}, oi);
OIs{2} = oiCompute(scene{2}, oi);

% Show irradiance (optical image) 
%vcAddObject(OIs{1}); oiWindow;
%vcAddObject(OIs{2}); oiWindow;

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