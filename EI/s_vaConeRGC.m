%% Starting with two oiSequences, create cone mosaic response
%

%%
ieInit

%% This script creates the two sequences, offset and aligned

% Check that they are about right.  The visualization needs work.

s_vaStimulus;
offset.visualize;
aligned.visualize;

%%
cMosaic = coneMosaic('os',osLinear);
% cMosaic = coneMosaic('os',osLinear,'size',[1 50]);
% cMosaic = coneMosaic('os',osBioPhys);

% These lines match the mosaic to the oiSequence
% The don't have to be matched, but maybe we should always do that.
fov = 0.6 * oiGet(offset.oiFixed,'fov');
cMosaic.setSizeToFOV(fov);

% Match the stimulus samples with the integration time
cMosaic.integrationTime = offset.timeAxis(2) - offset.timeAxis(1);

% Generate eye movements for the whole sequence
tSamples = length(offset.modulationFunction);
cMosaic.emGenSequence(tSamples);

%% Compute the absorptions, including eye movements

cMosaic.compute(offset);
cMosaic.name = 'offset';
cMosaic.window;

% cMosaic.compute(aligned);
% cMosaic.name = 'aligned';

%% Show the photocurrent, without or with noise

cMosaic.os.noiseFlag = true;
cMosaic.computeCurrent;
cMosaic.window;

%% Plot the impulse response

[~, meancurrent] = cMosaic.os.linearFilters(cMosaic);
cMosaic.plot('os current filters','meancurrent',meancurrent);

%% Testing
% 
% scene = sceneCreate('uniform ee');
% scene = sceneSet(scene,'fov',1);
% oi = oiCreate('human');
% oi = oiCompute(oi,scene);
% 
% %%
% cMosaic = coneMosaic;
% cMosaic.integrationTime = 0.001;
% cMosaic.setSizeToFOV(0.5);
% cMosaic.compute(oi);
% cMosaic.window;

%% Compute the bipolar response
 
bp = bipolar(cMosaic);
bp.set('sRFcenter',1);
bp.set('sRFsurround',1);
bp.compute(cMosaic);
 
% bp.plot('movie response')
 
%% Set RGC mosaic parameters
 
clear params innerRetinaSU
cellType = 'onParasol';
% cellType = 'offParasol';
params.name = 'macaque phys';
params.eyeSide = 'left'; 
params.eyeRadius = sqrt(sum(0.^2)); 
% params.fov = fov;
params.eyeAngle = 0; ntrials = 0;
 
% Create RGC object
innerRetinaSU = ir(bp, params);
innerRetinaSU.mosaicCreate('type',cellType,'model','GLM');
 
nTrials = 1; innerRetinaSU = irSet(innerRetinaSU,'numberTrials',nTrials);
 
%% Compute the inner retina response
 
innerRetinaSU = irCompute(innerRetinaSU, bp); 