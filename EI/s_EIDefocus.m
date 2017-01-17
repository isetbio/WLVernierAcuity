%% Impact of blurring the optics on CSF
%

% Show the dependence on spatial size of the cone mosaic for the computational
% observer.
nTrials = 1000;
nBasis  = 40;

% Integration time 
tStep   = 10;  % Adequate for absorptions (ms)  Should be in units of sec!

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.5;

% Original scene
sceneFOV = 1;

defocus = 1.5;

% Spatial scale to control visual angle of each display pixel The rule is 6/sc
% arc sec for a 0.35 deg scene. If you change the scene to 0.5 deg then 0.5/0.35
sc = 1*(sceneFOV/0.35);  

freqSamples = [1, 2, 10, 15, 20, 25];

s_EIParametersCSF;

% Special conditions
contrasts   = logspace(-2.5,0,5); 

%% Summarize
% 
% params.harmonic.contrast = 1;
% params.harmonic.freq = freqSamples(end);
% [~, harmonic,scenes,tseries] = csfStimuli(params);
% 
% % Show and ultimately print
% ieAddObject(scenes{2}); sceneWindow;
% ieAddObject(harmonic.oiModulated); oiWindow;
% degPerPixel = sceneGet(scenes{2},'degrees per sample');
% minPerPixel = degPerPixel * 60;
% secPerPixel = minPerPixel * 60;

%% Examine the oiSequence

% vcNewGraphWin;
% for ii=1:harmonic.length
%     oi = harmonic.frameAtIndex(ii);
%     imagesc(oiGet(oi,'rgb'));
%     title(sprintf('%d',ii)); pause(0.1);
% end

% Or examine the cone mosaic
%
% save blankFrame harmonic params
% 
% load('blankFrame')
% ieAddObject(harmonic.oiModulated); oiWindow;
% ieAddObject(harmonic.oiFixed); oiWindow;
% oi = harmonic.frameAtIndex(10); % ieAddObject(oi); oiWindow;
% vcNewGraphWin; imagesc(oiGet(oi,'rgb'));

%%
% cm = coneMosaic;
% cm.setSizeToFOV(params.cmFOV);
% cm.integrationTime = harmonic.timeStep;
% cm.noiseFlag = 'random';
% emPaths  = cm.emGenSequence(harmonic.length, 'nTrials', params.nTrials, ...
%     'em', params.em);
% cm.compute(harmonic);
% cm.window;

%% Main compute loop

% Could start the parallel pool
% c = gcp; if isempty(c), parpool('local'); end

tic;
PC = zeros(length(contrasts),length(freqSamples));
parfor pp=1:length(freqSamples)
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
lStrings = cell(1,length(freqSamples));
for pp=1:length(freqSamples)
    lStrings{pp} = sprintf('%.2f deg',freqSamples(pp));
end

h = vcNewGraphWin;
semilogx(contrasts,PC,'-o');
xlabel('Contrast'); ylabel('Percent correct')
grid on; l = legend(lStrings);
set(l,'FontSize',12)

%%
str = datestr(now,30);
fname = fullfile(wlvRootPath,'EI','figures',['csf-',str,'.mat']);
save(fname, 'PC','params', 'contrasts', 'freqSamples','scenes');

%%
