%% Initialize the base parameters for the EI paper
%  Typically, we run this before the segments
%

%%
nTrials = 300;

% Integration time and time step of the movie
tStep   = 30;  % Adequate for absorptions
% tStep   = 5;   % Useful for current

% Set the number of image bases
nBasis = 50;

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.25;

sc = 1;
% Scene field of view in degrees
sceneFOV      = 0.35;

%% Set basic parameters for the vernier stimulus
clear params;

v = vernierP;
v.gap = 0;
v.bgColor = 0.5;    % The dark background used in McKee and Westheimer

% For a scene fov of 0.35 and a size of 210,210, 1 pixel offset is 6 sec of
% arc.  Scaling them together preserves this 6 arc sec value. 
v.sceneSz = [210 210]*sc;
v.barWidth  = 10;
v.barLength = 200;

% Attach the vernier parameters
params.vernier = v;

%% These are oisequence and other parameters

params.tsamples  = (-200:tStep:200)*1e-3;    % In second M/W was 200 ms
params.timesd    = 100*1e-3;                 % In seconds                 
params.nTrials   = nTrials;
params.tStep     = tStep;
params.sc        = sc;
params.nBasis    = nBasis;
params.fov       = coneMosaicFOV;            % Cone mosaic field of view (deg)
params.distance  = 0.3;
params.em        = emCreate;
params.em.emFlag = [1 1 1]';

%%
