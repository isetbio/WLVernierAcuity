%% s_vaLayers.m
%
% Testing the new layer architecture.
%
% BW ISETBIO Team, 2016

%%
tic
coneMosaicData = fullfile(isetbioRootPath,'local','coneMosaicData.mat');

if exist(coneMosaicData,'file')
    load(coneMosaicData);
    disp('Loading saved cone mosaic data');
else
    
    % Show the dependence on the cone mosaic size for the computational
    % observer.
    nTrials = 10;
    nBasis  = 40;
    
    % Integration time
    tStep   = 10;         % Adequate for photocurrent (ms)
    
    % Cone mosaic field of view in degrees
    coneMosaicFOV = 0.35;
    
    % Original scene
    sceneFOV = 0.4;
    
    % Spatial scale to control visual angle of each display pixel The rule is
    % 6/sc arc sec for a 0.35 deg scene. If you change the scene to 0.5 deg
    % then 0.5/0.35
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
    % cMosaic.plot('impulse response'); cMosaic.plot('os impulse response');
    
    %% For aligned or offset
    
    disp('Computing cone mosaic eye movements');
    emPaths  = cMosaic.emGenSequence(tSamples, 'nTrials', nTrials, ...
        'em', params.em);
    
    % compute absorptions for aligned and offset
    disp('Computing cone mosaic current');
    [~,alignedC] = cMosaic.compute(aligned, 'currentFlag', true, ...
        'emPaths', emPaths);
    
    % Have a look cMosaic.window;
    
    
    % save(fullfile(isetbioRootPath,'local','coneMosaicData.mat'), 'cMosaic',
    % 'alignedC');
end


%% Create a set of bipolar cell types in the bipolar mosaic


% bp = bipolar(cMosaic,'cellType','onmidget');   % offdiffuse
% bp.set('sRFcenter',10); bp.set('sRFsurround',0);
%
% disp('Computing bipolar responses'); [~, bpNTrialsCenter,
% bpNTrialsSurround] = bp.compute(cMosaic,'coneTrials',alignedC);
%
% %

% Near fovea
patchEccentricity = 0;

% Make a cell of each type
cellType = {'ondiffuse','offdiffuse','onmidget','offmidget','onSBC'};
clear bpParams
bpParams.ecc = patchEccentricity;
bpParams.rectifyType = 1;

bpMosaic  = cell(1,length(cellType));
bpNTrials = cell(1,length(cellType));
for ii = 1:length(cellType)
    
    bpParams.cellType = cellType{ii};
    
    bpMosaic{ii} = bipolarMosaic(cMosaic, bpParams);
    bpMosaic{ii}.set('sRFcenter',1);
    bpMosaic{ii}.set('sRFsurround',0);
    
    [~, bpNTrialsCenterTemp, bpNTrialsSurroundTemp] = ...
        bpMosaic{ii}.compute(cMosaic,'coneTrials',alignedC);
    bpNTrials{ii} = bpNTrialsCenterTemp - bpNTrialsSurroundTemp;
    
end
%
% TODO:  bpMosaic object to make this more like innerRetina.mosaic{1}.
% After showing movie, the numbers on the mosaic axis are missing The units
% on the center size may not be correct. We need to allow changing the size
% of the center and surround on the bipolar. Maybe the size of the bipolar
% window got changed. Can we change the sizes of the bipolar receptive
% fields
%%
bpMosaic{1}.window;
%%
bpMosaic{2}.window;

%% Retinal ganlion cell model

clear rgcLayer irParams mosaicParams

% Choose a cell type
cellType = 'OFF Midget';  % 'offParasol'; 'onMidget' ...
irParams.name = 'macaque phys';
irParams.eyeSide = 'left';

% Create inner retina object
ecc = patchEccentricity;
irParams.eyeRadius = sqrt(sum(ecc.^2));
irParams.eyeAngle = 0; ntrials = 0;
rgcL = rgcLayer(bpMosaic, irParams);

% There are various parameters you could set.  We will write a script
% illustrating these later.  We need a description.
%
mosaicParams.centerNoise = 0;
mosaicParams.ellipseParams = [1 1 0];  % Principle, minor and theta
% mosaicParams.axisVariance = .1;
mosaicParams.type  = cellType;
mosaicParams.model = 'GLM';

diameters = [5 5 3 3 10];  % In microns.

cellType = {'on parasol','off parasol','on midget','off midget','smallbistratified'};
for ii = 1:length(cellType)
    mosaicParams.rfDiameter = diameters(ii);
    mosaicParams.type = cellType{ii};
    rgcL.mosaicCreate(mosaicParams);
end

nTrials = 1; rgcL.set('numberTrials',nTrials);
% innerRetina.mosaic{1}.get('rfDiameter')

%% Compute the inner retina response and visualize

% Number of trials refers to number of repeats of the same stimulus
disp('Computing rgc responses');
[rgcL, nTrialsSpikes] = rgcL.compute(bpMosaic,'bipolarTrials',bpNTrials);

%% Retinal ganglion cell layer window

% TODO:
%   Put up 'ieInWindowMessage() when the movie is playing Label the
%   distances on the x and y axes
%

% I wish we could have two windows up at the same time.  Read the code to
% see why it is always the same window.
rgcL.mosaic{1}.window;
rgcL.mosaic{3}.window;
rgcL.mosaic{5}.window;

%%

toc