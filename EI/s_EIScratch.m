%% Run s_vaAbsorptions

nTrials = 300;
% Integration time and time step of the movie
tStep   = 30;  % Adequate for absorptions
% tStep   = 5;   % Useful for current

% Set the number of image bases
nBasis = 15;

% Cone mosaic field of view in degrees
% The scene should always be less than 0.35 deg because that is the size in
% vaStimuli.
coneMosaicfov = 0.25;

% Sets the number of steps in the curve
barOffset = [4 6];

% Set basic parameters for the vernier stimulus
clear params;
v = vernierP;
v.gap = 0;
v.bgColor = 0.1;

% For a scene fov of 0.35 and this size, 1 pixel offset is 6 sec of arc
v.sceneSz = [210 210];
v.barWidth  = 10;
v.barLength = 100;
% Attach the vernier parameters
params.vernier = v;

% These are oisequence and other parameters
params.tsamples  = (-200:tStep:200)*1e-3;   % In second M/W was 200 ms
params.timesd  = 100*1e-3;                  % In seconds                 
params.nTrials = nTrials;
params.tStep   = tStep;
params.nBasis  = nBasis;
params.fov     = coneMosaicfov;             % Field of view of cone mosaic in deg
params.em      = emCreate;
params.em.emFlag = [1 1 1]';

% vaPCA(params);

% We need to visualize the imageBasis.  Here would be a good spot.
s_vaAbsorptions
% s_vaCurrent

%%
params = 

    barOffset: 6
     barWidth: 3
     tsamples: [1x22 double]
       timesd: 0.0400
      nTrials: 400
        tStep: 8

X = [0 2 4 6 ]
P = [51.25 88.12 97.50 100.00 ]
