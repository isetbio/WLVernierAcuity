%% Starting with two oiSequences, create cone mosaic response
%

%%
ieInit

%%  Create a brief flash for testing linearity

scene = sceneCreate('uniform ee');
oi = oiCreate;
oi = oiCompute(oi,scene);
levels = logspace(0,3,5);
OIs{1} = oiAdjustIlluminance(oi,0);
OIs{2} = oiAdjustIlluminance(oi,levels(1));
% for ii=1:2, ieAddObject(OIs{ii}); end; oiWindow;

integrationTime = 0.001;
sampleTimes = (1:150)*integrationTime;
nTimes = length(sampleTimes);

%%
cMosaic = coneMosaic('os',osLinear,'size',[1 3],'pattern',2 * ones(1,3));
cMosaic.integrationTime = integrationTime;
cMosaic.emGenSequence(nTimes);

%%
% The weights define some amount of the constant background and some amount
% of the line on the same constant background
weights = zeros(nTimes,1);
weights(50) = 1;
oiStep = oiSequence(OIs{1}, OIs{2}, ...
    sampleTimes, weights, ...
    'composition', 'add');
% oiStep.visualize('format','movie');

%%
cMosaic.compute(oiStep);

cMosaic.computeCurrent;

cMosaic.window;

cMosaic.absorptions



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
