%% Starting with two oiSequences, create cone mosaic response
%

%% This script creates the two sequences, oiSeqAligned and oiSeqOffset
s_vaStimulus;

%%
cMosaic = coneMosaic('os',osLinear);
% cMosaic = coneMosaic('os',osBioPhys);

cMosaic.setSizeToFOV(0.6 * imgFov);
cMosaic.integrationTime = 0.001;
cMosaic.emGenSequence(tSamples);

%% Compute the responses with eye movements
cMosaic.compute(oiSeqOffset);
cMosaic.name = 'offset';

%%
cMosaic.os.noiseFlag = false;
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
