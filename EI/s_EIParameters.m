%% Initialize the base parameters for the EI paper
%
% Typically, we run this before the segments.  Then we over-ride these values to
% produce different curves.
%
% We also set different parameters, such as the barOffset in pixels.
%
% See ...

%% These are oisequence and other parameters
clear params

params.tsamples  = (-200:tStep:200)*1e-3;   % In second M/W was 200 ms
params.timesd    = 100e-3;                  % +/- 1 sstd is 200 ms               
params.nTrials   = nTrials;
params.tStep     = tStep;
params.sc        = sc;
params.nBasis    = nBasis;
params.fov       = coneMosaicFOV;            % Cone mosaic field of view (deg)
params.distance  = 0.3;
params.em        = emCreate;
params.em.emFlag = [1 1 1]';


%% Set basic parameters for the vernier stimulus

v = vernierP;
v.gap     = 0;
v.bgColor = 0.5;    % The dark background used in McKee and Westheimer

% For a scene fov of 0.35 and a size of 210,210, 1 pixel offset is 6 sec of
% arc.  Scaling them together preserves this 6 arc sec value. 
v.sceneSz = [210 210]*params.sc;
v.barWidth  = 10;
v.barLength = 200;

% Attach the vernier parameters
params.vernier = v;

%%
