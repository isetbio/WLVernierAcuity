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
params.cmFOV     = coneMosaicFOV;           % Cone mosaic field of view (deg)
params.sceneFOV  = sceneFOV;                % Scene field of view (deg)
params.distance  = 0.3;
params.em        = emCreate;
params.em.emFlag = [1 1 1]';


%% Set basic parameters for the harmonic stimulus

h = harmonicP;
h.row = 210*params.sc;
h.col = 210*params.sc;
h.GaborFlag = 0.2;

% Attach the vernier parameters
params.harmonic = h;

%%
