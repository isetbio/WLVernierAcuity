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
tStep   = 10;         % Adequate for photocurrent (ms)

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
%%
% Got the params from elsewhere ...
%
%  load ...
params.vernier.offset = 0;

[aligned, offset, ~, ~] = vaStimuli(params);

%  Compute absorptions for multiple trials
tSamples = aligned.length;
cMosaic = coneMosaic;

% Sometimes we set the mosaic size to 15 minutes (.25 deg) because that is
% the spatial pooling size found by Westheimer and McKee
cMosaic.setSizeToFOV(params.cmFOV);

% Not sure why these have to match, but there is a bug if they don't.
cMosaic.integrationTime = aligned.timeStep;

% For aligned or offset
cMosaic.noiseFlag = 'random';
emPaths  = cMosaic.emGenSequence(tSamples, 'nTrials', nTrials, ...
    'em', params.em);

% compute absorptions for aligned and offset
[~,alignedC] = cMosaic.compute(aligned, 'currentFlag', true, ...
    'emPaths', emPaths);

% Have a look
% cMosaic.window;

%%
bp = bipolar(cMosaic);
bp.set('sRFcenter',10);
bp.set('sRFsurround',0);

%
[~, bpNTrialsCenter, bpNTrialsSurround] = bp.compute(cMosaic,'coneTrials',alignedC);

% % Should make a loop option for the movie window;
% % Should have the movie window force a new window with a stop button.
% vcNewGraphWin;
% ieMovie(squeeze(bpNTrialsCenter(1,:,:,:)));
bp.window;
%%
clear innerRetina

% Choose a cell type
cellType = 'onParasol';
% cellType = 'offParasol';
rgcparams.name = 'macaque phys';
rgcparams.eyeSide = 'left';

ecc = 0;
rgcparams.eyeRadius = sqrt(sum(ecc.^2)); 
rgcparams.eyeAngle = 0; ntrials = 0;
 
% Create RGC object
% TODO: Do the size/position randomness with a function like emGenerate();
% That would take in (rows,cols) and noise parameters to generate the set
% of ellipse parameters.
innerRetina = ir(bp, rgcparams);
rgcParams.type = cellType;
rgcParams.centerNoise = .2;
rgcParams.model = 'GLM';
rgcParams.ellipseParams = []; %[1 1 0];  % Principle, minor and theta
innerRetina.mosaicCreate(rgcParams);
innerRetina.mosaic{1}.window;

% innerRetina.mosaicCreate('type',cellType,'model','GLM','centerNoise',centerNoise);

% Number of trials refers to number of repeats of the same stimulus
nTrials = 1; innerRetina.set('numberTrials',nTrials);
 
%% Compute the inner retina response

% 
[innerRetina, nTrialsSpikes] = innerRetina.compute(bp,'bipolarTrials',bpNTrialsCenter - bpNTrialsSurround); 
lastTime = innerRetina.mosaic{1}.get('last spike time');
 
%% Make the PSTH movie
% innerRetina.mosaic{1}.set('dt',1);
% psth = innerRetina.mosaic{1}.get('psth');
% allCells = mean(mean(psth,1),2);
% allCells = squeeze(allCells);
% vcNewGraphWin; plot(allCells)

% Let's make this work.
% vcNewGraphWin;
% innerRetina.mosaic{1}.plot('psth');
% or ...
% innerRetina.plot('psth','mosaic',val);

clear movieparams 
movieparams.FrameRate = 5; movieparams.step = 2; movieparams.show = true;
 
% % View movie of RGC linear response
vcNewGraphWin; ieMovie(innerRetina.mosaic{1}.responseLinear);
 
% View movie of PSTH for mosaic
steadyStateFrame = 50; % Get rid of transient spiking
psth = innerRetina.mosaic{1}.get('psth');
vcNewGraphWin; movieparams.step = 10;ieMovie(psth,movieparams);

%% See the data

% Could be - 
%   innerRetina.window{mosaicNumber);
innerRetina.mosaic{1}.window;

%%