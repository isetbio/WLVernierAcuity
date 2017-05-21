%% s_vaRGC.m
%
% First experiment with RGC classification
%
% Deprecate the bottom
%
% HJ, ISETBIO TEAm, 2016

%%
% Show the dependence on the cone mosaic size for the computational
% observer.
nTrials = 2;
nBasis  = 40;

% Integration time 
tStep   = 5;         % Adequate for photocurrent (ms)

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.35;

% Original scene
sceneFOV = 0.4;

% Spatial scale to control visual angle of each display pixel The rule is 6/sc
% arc sec for a 0.35 deg scene. If you change the scene to 0.5 deg then 0.5/0.35
sc = 3*(sceneFOV/0.35);   % If you do not multiply by a scalar, offset is 6 arc sec

s_EIParameters;

% Make the bar length a little less than the scene size
params.vernier.barLength = params.vernier.sceneSz(1)-1;
params.tsamples  = (-200:tStep:400)*1e-3;
%% Read the stimulus, which might have been saved

params.vernier.offset = 0;
[aligned, offset, ~, ~] = vaStimuli(params);

%%  Compute absorptions for multiple trials
tSamples = aligned.length;
cMosaic = coneMosaic;

% Sometimes we set the mosaic size to 15 minutes (.25 deg) because that is
% the spatial pooling size found by Westheimer and McKee
cMosaic.setSizeToFOV(params.cmFOV);

% Not sure why these have to match, but there is a bug if they don't.
cMosaic.integrationTime = aligned.timeStep;
cMosaic.noiseFlag = 'random';
% cMosaic.plot('impulse response');
% cMosaic.plot('os impulse response');

%% For aligned or offset

disp('Computing cone mosaic eye movements');
emPaths  = cMosaic.emGenSequence(tSamples, 'nTrials', nTrials, ...
    'em', params.em);

% compute absorptions for aligned and offset
disp('Computing cone mosaic current');
[~,alignedC] = cMosaic.compute(aligned, 'currentFlag', true, ...
    'emPaths', emPaths);

% Have a look
% cMosaic.window;

%%
bp = bipolar(cMosaic,'cellType','onmidget');   % offdiffuse
bp.set('sRFcenter',10);
bp.set('sRFsurround',0);

disp('Computing bipolar responses');
[~, bpNTrialsCenter, bpNTrialsSurround] = bp.compute(cMosaic,'coneTrials',alignedC);

% bp.window;
%% Retinal ganlion cell model

clear innerRetina irParams mosaicParams innerRetina

% Choose a cell type
cellType = 'OFF Midget';  % 'offParasol'; 'onMidget' ...
irParams.name = 'macaque phys';
irParams.eyeSide = 'left';

% Create inner retina object
ecc = 0;
irParams.eyeRadius = sqrt(sum(ecc.^2)); 
irParams.eyeAngle = 0; ntrials = 0;
irParams.eyeRadius = 0;
innerRetina = ir(bp, irParams);

% There are various parameters you could set.  We will write a script
% illustrating these later.  We need a description.
%
mosaicParams.centerNoise = 0.2;
% mosaicParams.ellipseParams = [1 .8 0];  % Principle, minor and theta
% mosaicParams.axisVariance = .1;
mosaicParams.type  = cellType;
mosaicParams.model = 'GLM';
innerRetina.mosaicCreate(mosaicParams);

% innerRetina.mosaic{1}.set('rfDiameter',10);
innerRetina.mosaic{1}.rgcInitSpace(innerRetina,cellType);

nTrials = 1; innerRetina.set('numberTrials',nTrials);
innerRetina.mosaic{1}.get('rfDiameter')

%% Compute the inner retina response and visualize

% Number of trials refers to number of repeats of the same stimulus
disp('Computing rgc responses');
[innerRetina, nTrialsSpikes] = innerRetina.compute(bp,'bipolarTrials',bpNTrialsCenter - bpNTrialsSurround); 
 
% Could become - innerRetina.window{mosaicNumber);
innerRetina.mosaic{1}.window;

%% 