%% t_VernierCones
%
% Illustrate the parameters of human cones
%
%  HJ/BW, ISETBIO TEAM, 2015

%% Initialize a new ISETBIO session
ieInit;

%% Create a verneir scene
%  See t_VernierScene for options in vernier scene create. Here, we use a
%  basic version, composing of two misaligned lines
p.bgColor = 0.5;
p.offset  = 1;
scene = sceneCreate('vernier', 'display', p);
vcAddObject(scene); sceneWindow;

%% Create standard human optics
%  See t_VernierOptics for options in creating human optics. Here, we
%  create a standard one with wavefront measurements
oi = oiCreate('wvf human');

% compute irradiance
oi = oiCompute(oi, scene);
vcAddObject(oi); oiWindow;

%% Create standard human cone mosaic
sensor = sensorCreate('human');

% adjust cone mosaic size
sensor = sensorSetSizeToFOV(sensor, sceneGet(scene, 'fov'), scene, oi);

% compute cone absorptions
sensor = sensorCompute(sensor,oi);

% have a look
vcAddObject(sensor); sensorWindow('scale',1);

%% Create customized human cone mosaic
%  create cone with customized K,L,M,S density
cone = coneCreate('human', 'spatial density', [0 0.3 0.6 0.1]);

%  adjust parameters
cone = coneSet(cone, 'peak efficiency', [0.6 0.6 0.6]);

%  create coen mosaic
sensor = sensorCreate('human', cone);

%  adjust cone size
sensor = sensorSet(sensor, 'pixel pd width', 3e-6);
sensor = sensorSet(sensor, 'pixel pd height', 3e-6);

%  adjust cone temporal integration time
sensor = sensorSet(sensor, 'exp time', 0.04);

%  compute cone absorptions with customized sensor
sensor = sensorCompute(sensor, oi);
vcAddObject(sensor); sensorWindow('scale', 1);

%% Eye-movement simulation
% Number of exposures samples
nFrames = 50;

expTime = sensorGet(sensor, 'exp time');   % Usually about 50 ms
emDuration = 0.001;
emPerExposure = expTime / emDuration;
sensor = sensorSet(sensor, 'exp time', emDuration);

% Generate eyemovement
p.nSamples = nFrames * emPerExposure;
p.emFlag   = [1 1 1];       % Include tremor drift and saccade
sensor = eyemoveInit(sensor, p);

% The coneAbsorptions function is an interface to sensorCompute. Notice
% that when we make an eye movement video we call coneAbsorptions, not
% sensorCompute.
cones = coneAbsorptions(sensor, oi);

% Show the eye movement positions
ePos = sensorGet(cones,'sensor positions');
vcNewGraphWin;
plot(ePos(:,1),ePos(:,2),'-o')
xlabel('Cone position');
ylabel('Cone Position')
grid on

%% END