%% s_EIScratch
% Sets up parameters and runs s_vaAbsorptions

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
nTrials = 300;

% Integration time 
tStep   = 30;  % Adequate for absorptions (ms)

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.25;

% Spatial scale to control visual angle of each display pixel
sc = 1;

s_EIParameters;
barOffset = [0 1  3  5 ];
vals = [10 40 120 180 240];
PC = zeros(length(barOffset),length(vals));

for pp=1:length(vals)
    params.vernier.barLength = vals(pp);
    s_vaAbsorptions;
    PC(:,pp) = P(:);
end
% mesh(PC)

% Offset per sample on the display
offsetDeg = sceneGet(scenes{1},'degrees per sample');
offsetSec = offsetDeg*3600;

% Legend
lStrings = cell(1,length(vals));
for pp=1:length(vals)
    lStrings{pp} = sprintf('%.1f arc min',offsetSec*vals(pp)/60);
end

h = vcNewGraphWin;
plot(offsetSec*barOffset,PC);
xlabel('Offset arc sec'); ylabel('Percent correct')
grid on; l = legend(lStrings);
set(l,'FontSize',12)

fname = fullfile(wlvRootPath,'EI','figures','spatialBarLength.mat');
save(fname, 'PC','params', 'barOffset', 'vals');
 
%%  Sweep out different durations
%

% Total of 1 sec duration
params.tsamples  = (-500:tStep:500)*1e-3;   % In second M/W was 200 ms
barOffset = 1;
sd = [50 100 200 400 600]*1e-3;
PC = zeros(1,length(sd));

for pp=1:length(sd)
    params.timesd  = sd(pp);                  % In seconds                 
    s_vaAbsorptions;
    PC(:,pp) = P(:);
end

% Offset per sample on the display
offsetDeg = sceneGet(scenes{1},'degrees per sample');
offsetSec = offsetDeg*3600*barOffset;

vcNewGraphWin;
plot(sd,squeeze(PC));
xlabel('Duration (ms)'); ylabel('Percent correct')
grid on; l = legend(sprintf('Offset: %d arc sec',offsetSec));
set(l,'FontSize',12);

fname = fullfile(wlvRootPath,'EI','figures','temporalPooling.mat');
save(fname, 'PC', 'params', 'barOffset', 'sd');

%%  Fix the eye movements to 0.
% Purpose:  Show that eye movements have a big impact
%

nTrials = 300;

% Integration time 
tStep   = 30;  % Adequate for absorptions (ms)

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.25;

% Spatial scale to control visual angle of each display pixel
sc = 1;

s_EIParameters;

params.vernier.barLength = 120;

% Maybe compare the prob. correct at 6 sec when there are no eye movements to
% the standard eye movement parameter
barOffset = [0 1 3];
PC = zeros(numel(barOffset),2);
params.em.emFlag = [0 0 0]';
s_vaAbsorptions;
PC(:,1) = P(:);

params.em.emFlag = [1 1 1]';
s_vaAbsorptions;
PC(:,2) = P(:);

% Offset per sample on the display
offsetDeg = sceneGet(scenes{1},'degrees per sample');
offsetMin = offsetDeg*60;

% Legend
lStrings = cell({'No em','Default em'});

vcNewGraphWin;
plot(offsetMin*barOffset*60,PC);
xlabel('Offset arc sec'); ylabel('Percent correct')
grid on; l = legend(lStrings);
set(l,'FontSize',12)

fname = fullfile(wlvRootPath,'EI','figures','eyeMovements.mat');
save(fname, 'PC', 'params', 'barOffset');


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
