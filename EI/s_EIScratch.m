%% s_EIScratch
% Sets up parameters and runs s_vaAbsorptions

nTrials = 300;

% Integration time and time step of the movie
tStep   = 30;  % Adequate for absorptions
% tStep   = 5;   % Useful for current

% Set the number of image bases
nBasis = 40;

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.25;

sc = 1;
% Scene field of view in degrees
sceneFOV      = 0.35;

% Sets the number of steps in the curve
% barOffset = [0 2 4 6];

% Set basic parameters for the vernier stimulus
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

% These are oisequence and other parameters
params.tsamples  = (-200:tStep:200)*1e-3;   % In second M/W was 200 ms
params.timesd  = 100*1e-3;                  % In seconds                 
params.nTrials = nTrials;
params.tStep   = tStep;
params.sc      = sc;
params.nBasis  = nBasis;
params.fov     = coneMosaicFOV;             % Field of view of cone mosaic in deg
params.distance  = 0.3;
params.em        = emCreate;
params.em.emFlag = [1 1 1]';

% vaPCA(params);

% We need to visualize the imageBasis.  Here would be a good spot.
%
% s_vaAbsorptions
% s_vaAbsorptionsHex

% s_vaCurrent

%% Loop on bar length

% Show the dependence on bar length for the computational observer.
% Use the match on bar length as an indicator of the spatial summation region of
% the human eye

barOffset = [0 1 2 3 4 5 6];
vals = [10 40 160];    % 5 arc sec per pixel
PC = zeros(length(barOffset),length(vals));

for pp=1:length(vals)
    params.vernier.barLength = vals(pp);
    s_vaAbsorptions;
    PC(:,pp) = P(:);
end

% Offset per sample on the display
offsetDeg = sceneGet(scenes{1},'degrees per sample');
offsetSec = offsetDeg*3600;

% Legend
lStrings = cell(1,length(vals));
for pp=1:length(vals)
    lStrings{pp} = sprintf('%.1f arc sec',offsetSec*vals(pp));
end

h = vcNewGraphWin;
plot(offsetSec*barOffset,PC);
xlabel('Offset arc sec'); ylabel('Percent correct')
grid on; l = legend(lStrings);
set(l,'FontSize',12)

fname = fullfile(wlvRootPath,'EI','figures','spatialBarLength.mat');
save(fname, 'params', 'barOffset', 'vals');
 
%%  Sweep out different durations
%
% Run with sc = 1

params.tsamples  = (-500:tStep:500)*1e-3;   % In second M/W was 200 ms

barOffset = 2;
sd = [50 100 200 400]*1e-3;
PC = zeros(1,length(sd));
for pp=1:length(sd)
    params.timesd  = sd(pp);                  % In seconds                 
    s_vaAbsorptions;
    PC(:,pp) = P(:);
end

% Offset per sample on the display
offsetDeg = sceneGet(scenes{1},'degrees per sample');
offsetSec = offsetDeg*3600;

vcNewGraphWin;
plot(sd,squeeze(PC));
xlabel('Duration (ms)'); ylabel('Percent correct')
grid on; l = legend(sprintf('Offset: %d arc sec',offsetSec));
set(l,'FontSize',12);

params

%%  Fix the eye movements to 0.
% Purpose:  Show that eye movements have a big impact
%

% Maybe compare the prob. correct at 6 sec when there are no eye movements to
% the standard eye movement parameter
params.em.emFlag = [0 0 0]';
barOffset = [0 1 2 4 6];

s_vaAbsorptions;

% Offset per sample on the display
offsetDeg = sceneGet(scenes{1},'degrees per sample');
offsetSec = offsetDeg*3600;

% Legend
lStrings = cell(1,length(barOffset));
for pp=1:length(barOffset)
    lStrings{pp} = sprintf('%.1f arc sec',offsetSec*barOffset(pp));
end

vcNewGraphWin;
plot(offsetSec*barOffset,P);
xlabel('Offset arc sec'); ylabel('Percent correct')
grid on; l = legend(lStrings);
set(l,'FontSize',12)

%% Vary the spatial blurring by the optics 
%
% Purpose: Show that blurring impacts the ability to resolve a 30 c/deg target, but does
% not impact the vernier resolution




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
