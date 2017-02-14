%% s_vaRGC.m
%
% First experiment with RGC classification
%
% Deprecate the bottom
%
% HJ, ISETBIO TEAm, 2016

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
% cMosaic.window;

bp = bipolar(cMosaic);
bp.set('sRFcenter',10);
bp.set('sRFsurround',0);
[~, bpNTrialsCenter, bpNTrialsSurround] = bp.compute(cMosaic,'nTrialsInput',alignedC);

%%
clear innerRetinaSU
cellType = 'onParasol';
% cellType = 'offParasol';
rgcparams.name = 'macaque phys';
rgcparams.eyeSide = 'left';

ecc = 0;
rgcparams.eyeRadius = sqrt(sum(ecc.^2)); 
rgcparams.eyeAngle = 0; ntrials = 0;
 
% Create RGC object
innerRetinaSU = ir(bp, rgcparams);
innerRetinaSU.mosaicCreate('type',cellType,'model','GLM');
 
nTrials = 1; innerRetinaSU = irSet(innerRetinaSU,'numberTrials',nTrials);
 
%% Compute the inner retina response
 
[innerRetinaSU, nTrialsSpikes] = irCompute(innerRetinaSU, bp,'nTrialsInput',bpNTrialsCenter-bpNTrialsSurround); 
lastTime = innerRetinaSU.mosaic{1}.get('last spike time');
 
%% Make the PSTH movie
innerRetinaSU.mosaic{1}.set('dt',1);
psth = innerRetinaSU.mosaic{1}.get('psth');
 
clear movieparams 
movieparams.FrameRate = 5; movieparams.step = 2; movieparams.show = true;
 
% % View movie of RGC linear response
vcNewGraphWin; ieMovie(innerRetinaSU.mosaic{1}.responseLinear);
 
% View movie of PSTH for mosaic
steadyStateFrame = 50; % Get rid of transient spiking
vcNewGraphWin; ieMovie(psth,movieparams);
% [~,offsetC] = cMosaic.compute(offset, 'currentFlag', true, ...
%     'emPaths', emPaths);

