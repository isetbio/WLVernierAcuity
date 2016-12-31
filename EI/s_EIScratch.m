%% s_EIScratch
% Sets up parameters and runs s_vaAbsorptions

nTrials = 500;

% Integration time and time step of the movie
tStep   = 20;  % Adequate for absorptions
% tStep   = 5;   % Useful for current

% Set the number of image bases
nBasis = 40;

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.25;

sc = 2;
% Scene field of view in degrees
sceneFOV      = 0.35;
% Sets the number of steps in the curve
barOffset = [0 2 4 6];

% Set basic parameters for the vernier stimulus
clear params;
v = vernierP;
v.gap = 0;
v.bgColor = 0.02;    % The dark background used in McKee and Westheimer

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
params.nBasis  = nBasis;
params.fov     = coneMosaicFOV;             % Field of view of cone mosaic in deg
params.distance = 0.3;
params.em      = emCreate;
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
vals = [25 50 100 150 200];
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

vcNewGraphWin;
plot(offsetSec*barOffset,PC);
xlabel('Offset arc sec'); ylabel('Percent correct')
grid on; l = legend(lStrings);
set(l,'FontSize',12)

%%  Sweep out different durations

params.tsamples  = (-500:tStep:500)*1e-3;   % In second M/W was 200 ms

barOffset = 2;
sd = [20 50 100 200 400]*1e-3;
PC = zeros(1,length(sd));
for pp=1:length(sd)
    params.timesd  = sd(pp);                  % In seconds                 
    s_vaAbsorptions;
    PC(:,pp) = P(:);
end

% Offset per sample on the display
offsetDeg = sceneGet(scenes{1},'degrees per sample');
offsetSec = offsetDeg*3600;

% Legend
lStrings = cell(1,length(vals));
for pp=1:length(vals)
    lStrings{pp} = sprintf('%.1f arc sec',sd(pp));
end

vcNewGraphWin;
plot(sd,squeeze(PC));
xlabel('Duration (ms)'); ylabel('Percent correct')
grid on; l = legend(lStrings);
set(l,'FontSize',12)

%%  Fix the eye movements to 0.

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
