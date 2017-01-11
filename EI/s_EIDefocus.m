%% Impact of blurring the optics on CSF and Vernier
%

% Show the dependence on spatial size of the cone mosaic for the computational
% observer.
nTrials = 100;
nBasis  = 40;

% Integration time 
tStep   = 10;  % Adequate for absorptions (ms)  Should be in units of sec!

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.5;

% Original scene
sceneFOV = 1;

defocus = 0.5;

% Spatial scale to control visual angle of each display pixel The rule is 6/sc
% arc sec for a 0.35 deg scene. If you change the scene to 0.5 deg then 0.5/0.35
sc = 1*(sceneFOV/0.35);  

freqSamples = 1;  % CPD

s_EIParametersCSF;

%% Summarize

params.harmonic.freq = freqSamples(1);
params.harmonic.contrast = 1;
[~, harmonic,scenes,tseries] = csfStimuli(params);

% Show and ultimately print
ieAddObject(scenes{2}); sceneWindow;
ieAddObject(offset.oiModulated); oiWindow;
degPerPixel = sceneGet(scenes{2},'degrees per sample');
minPerPixel = degPerPixel * 60;
secPerPixel = minPerPixel * 60;

%% Set up for the CSF calculation

PC = zeros(length(contrasts),length(freqSamples));

%% 
% c = gcp; if isempty(c), parpool('local'); end

contrasts = 0.5;

tic;
for pp=1:length(freqSamples)
    fprintf('Starting %d of %d ...\n',pp,length(freqSamples));
    thisParam = params;
    thisParam.harmonic.freq = freqSamples(pp);
    P = csfAbsorptions(contrasts,thisParam);
    PC(:,pp) = P(:);
    fprintf('Finished %d\n',pp);
end
toc


%% Make summary graphs

% Legend
lStrings = cell(1,length(cmFOV));
for pp=1:length(cmFOV)
    lStrings{pp} = sprintf('%.2f deg',cmFOV(pp));
end

h = vcNewGraphWin;
plot(secPerPixel*barOffset,PC);
xlabel('Offset arc sec'); ylabel('Percent correct')
grid on; l = legend(lStrings);
set(l,'FontSize',12)

%%
str = datestr(now,30);
fname = fullfile(wlvRootPath,'EI','figures',['mosaicSize-',str,'.mat']);
save(fname, 'PC','params', 'contrasts', 'cmFOV','scenes');

%%
